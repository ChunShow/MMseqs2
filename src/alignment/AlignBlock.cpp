#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
// #include "CovSeqidQscPercMinDiag.lib.h"
// #include "CovSeqidQscPercMinDiagTargetCov.lib.h"
#include "QueryMatcher.h"
#include "IndexReader.h"
#include "FastSort.h"
#include "BlockAligner.h"
#include "Alignment.h"
#include "AlignmentSymmetry.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <set>  

#define MAX_SIZE 4096
#define MIN_SIZE 32

// Result buffer structure
struct ClusterResult {
    std::vector<unsigned int> memberIds;
    bool isValid;
    
    ClusterResult() : isValid(false) {}
};


static float parsePrecisionLib(const std::string &scoreFile, double targetSeqid, double targetCov, double targetPrecision) {
    std::stringstream in(scoreFile);
    std::string line;
    // find closest lower seq. id in a grid of size 5
    int intTargetSeqid = static_cast<int>((targetSeqid + 0.0001) * 100);
    int seqIdRest = (intTargetSeqid % 5);
    targetSeqid = static_cast<float>(intTargetSeqid - seqIdRest) / 100;
    // find closest lower cov. id in a grid of size 10
    targetCov = static_cast<float>(static_cast<int>((targetCov + 0.0001) * 10)) / 10;
    while (std::getline(in, line)) {
        std::vector<std::string> values = Util::split(line, " ");
        float cov = strtod(values[0].c_str(), NULL);
        float seqid = strtod(values[1].c_str(), NULL);
        float scorePerCol = strtod(values[2].c_str(), NULL);
        float precision = strtod(values[3].c_str(), NULL);
        if (MathUtil::AreSame(cov, targetCov) && MathUtil::AreSame(seqid, targetSeqid) && precision >= targetPrecision) {
            return scorePerCol;
        }
    }
    Debug(Debug::WARNING) << "Can not find any score per column for coverage "
                          << targetCov << " and sequence identity " << targetSeqid << ". No hit will be filtered.\n";

    return 0;
}

static void writeData(DBWriter *dbw, const std::pair<unsigned int, unsigned int> * ret, size_t dbSize) {
    std::string resultStr;
    resultStr.reserve(1024*1024*1024);
    char buffer[32];
    unsigned int prevRepresentativeKey = UINT_MAX;
    for(size_t i = 0; i < dbSize; i++){
        unsigned int currRepresentativeKey = ret[i].first;
        // write query key first
        if(prevRepresentativeKey != currRepresentativeKey) {
            if(prevRepresentativeKey != UINT_MAX){ // skip first
                dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
            }
            resultStr.clear();
            char *outpos = Itoa::u32toa_sse2(currRepresentativeKey, buffer);
            resultStr.append(buffer, (outpos - buffer - 1));
            resultStr.push_back('\n');
        }
        unsigned int memberKey = ret[i].second;
        if(memberKey != currRepresentativeKey){
            char * outpos = Itoa::u32toa_sse2(memberKey, buffer);
            resultStr.append(buffer, (outpos - buffer - 1) );
            resultStr.push_back('\n');
        }

        prevRepresentativeKey = currRepresentativeKey;
    }
    if(prevRepresentativeKey != UINT_MAX){
        dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
    }
}

int doAlign2clust(Parameters &par,
    DBWriter &resultWriter,
    DBReader<unsigned int> &alnDbr,
    const size_t dbFrom, const size_t dbSize) {
    DBReader<unsigned int> * seqDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);

    BaseMatrix *subMat;
    subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    int gapOpen = par.gapOpen.values.aminoacid();
    int gapExtend = par.gapExtend.values.aminoacid();
    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);

    float scorePerColThr = 0.0;
    std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
                                    
    scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    std::cout << "Score per column threshold for filtering: " << scorePerColThr << std::endl;
    EvalueComputation evaluer(seqDbr->getAminoAcidDBSize(), subMat);
    int32_t x_drop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());
    
    // Setup linear data
    unsigned int *assignedcluster = new(std::nothrow) unsigned int[dbSize];
    Util::checkAllocation(assignedcluster, "Can not allocate assignedcluster memory in Align2Clust");
    std::fill_n(assignedcluster, dbSize, UINT_MAX);


    // Determine optimal chunk size based on dataset size
    size_t totalItems = alnDbr.getSize();
    size_t numThreads = static_cast<size_t>(par.threads);

    // 총 항목의 일정 비율을 청크로 설정 (예: 1-5%)
    double chunkRatio = 0.02; // 2%
    size_t chunkSize = std::max(
        static_cast<size_t>(totalItems * chunkRatio),
        numThreads * 10
    );
    
    std::cout << "Using chunk size: " << chunkSize << " for " << totalItems << " items with " << par.threads << " threads" << std::endl;
    
    // Allocate result buffer for one chunk
    std::vector<ClusterResult> resultBuffer(chunkSize);

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
        Sequence query(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        Sequence target(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        BlockAligner blockaligner(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &fastMatrix, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
        char buffer[1024 + 32768*4];
        std::vector<std::pair<unsigned int, unsigned short>> dbKeysNdiagonal;
        dbKeysNdiagonal.reserve(1000); // Pre-allocate

        for (size_t chunkStart = 0; chunkStart < alnDbr.getSize(); chunkStart += chunkSize) {
            size_t chunkEnd = std::min(chunkStart + chunkSize, alnDbr.getSize());
            size_t availChunkSize = chunkEnd - chunkStart;
            
            // Parallel processing with dynamic scheduling and nowait
#pragma omp for schedule(dynamic, 1) nowait
            for (size_t i = chunkStart; i < chunkEnd; i++) {
                size_t bufferIdx = i - chunkStart;
                
                // Initialize result buffer for this item
                resultBuffer[bufferIdx].memberIds.clear();
                resultBuffer[bufferIdx].isValid = false;
                
                dbKeysNdiagonal.clear();
                
                const unsigned int queryKey = seqDbr->getDbKey(i);
                const size_t alnId = alnDbr.getId(queryKey);
                char *data = alnDbr.getData(alnId, thread_idx);
                size_t representative = seqDbr->getId(queryKey);
                size_t queryId = representative;
                
                // Early skip if already assigned
                if (assignedcluster[representative] != UINT_MAX) {
                    resultBuffer[bufferIdx].memberIds.push_back(UINT_MAX); // Mark as already assigned
                    continue;
                }
                
                // Mark as valid and add self-reference
                resultBuffer[bufferIdx].isValid = true;
                
                // Get query sequence
                char* querySequence = seqDbr->getData(queryId, thread_idx);
                size_t queryLen = seqDbr->getSeqLen(queryId);
                query.mapSequence(queryId, queryKey, querySequence, queryLen);
                blockaligner.initQuery(&query);
                matcher.initQuery(&query);
                
                // Parse prefilter hits
                while (*data != '\0') {
                    hit_t result = QueryMatcher::parsePrefilterHit(data);
                    const size_t targetId = seqDbr->getId(result.seqId);
                    if (assignedcluster[targetId] == UINT_MAX) {
                        dbKeysNdiagonal.push_back(std::make_pair(result.seqId, result.diagonal)); // query + 모든 targets
                    }
                    data = Util::skipLine(data);
                }

                // Process each target sequence
                for (size_t id = 0; id < dbKeysNdiagonal.size(); id++) {
                    const unsigned int targetKey = dbKeysNdiagonal[id].first;
                    const unsigned short diagonal = dbKeysNdiagonal[id].second;
                    const size_t targetId = seqDbr->getId(targetKey);
                    
                    // Skip identity matches
                    const bool isIdentity = (queryKey == targetKey);
                    if (isIdentity) {
                        resultBuffer[bufferIdx].memberIds.push_back(queryId);
                        continue;
                    }
                    
                    char* targetSequence = seqDbr->getData(targetId, thread_idx);
                    size_t targetLen = seqDbr->getSeqLen(targetId);
                    target.mapSequence(targetId, targetKey, targetSequence, targetLen);

                    // 1. Run ungapped alignment
                    BlockAligner::UngappedAln_res ungapped_aln = blockaligner.ungappedAlign(&target, diagonal); 
                    
                    // Check ungapped criteria
                    bool hasEvalue = (ungapped_aln.eval <= par.evalThr);
                    bool hasAlnLen = (ungapped_aln.alnLen >= par.alnLenThr);
                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, ungapped_aln.qcov, ungapped_aln.tcov);
                    float seqId = 0;
                    
                    if (hasEvalue) {    
                        int idCnt = 0;
                        for (int q = ungapped_aln.qStart; q <= ungapped_aln.qEnd; q++) {
                            char qLetter = querySequence[q] & static_cast<unsigned char>(~0x20);
                            char tLetter = targetSequence[ungapped_aln.tStart + (q - ungapped_aln.qStart)] & static_cast<unsigned char>(~0x20);
                            idCnt += (qLetter == tLetter) ? 1 : 0;
                        }
                        seqId = Util::computeSeqId(par.seqIdMode, idCnt, query.L, target.L, ungapped_aln.alnLen);
                    }
                    
                    char *end = Itoa::i32toa_sse2(ungapped_aln.alnLen, buffer);
                    size_t len = end - buffer;
                    std::string backtrace = "";
                    if (par.addBacktrace) {
                        backtrace = std::string(buffer, len - 1);
                        backtrace.push_back('M');
                    }
                    
                    if (isIdentity) {
                        ungapped_aln.qcov = 1.0f;
                        ungapped_aln.tcov = 1.0f;
                        seqId = 1.0f;
                    }
                    
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    
                    // If ungapped alignment passes all criteria
                    if (isIdentity || (hasAlnLen && hasCov && hasSeqId && hasEvalue)) {
                        Matcher::result_t result;
                        result = Matcher::result_t(targetKey, ungapped_aln.bitScore, ungapped_aln.qcov, ungapped_aln.tcov, seqId, ungapped_aln.eval, ungapped_aln.alnLen,
                                                   ungapped_aln.qStart, ungapped_aln.qEnd, query.L, ungapped_aln.tStart, ungapped_aln.tEnd, target.L, backtrace);
                        resultBuffer[bufferIdx].memberIds.push_back(targetId);
                        continue;
                    }
                    
                    // 2. Check if gapped alignment is needed
                    float currScorePerCol = static_cast<float>(ungapped_aln.score) / static_cast<float>(ungapped_aln.diagonalLen);
                    if (currScorePerCol < scorePerColThr) {
                        continue;
                    }
                    
                    // 3. Run gapped alignment
                    int alnLen = ungapped_aln.alnLen;
                    int qStartPos = ungapped_aln.qStart;
                    int tStartPos = ungapped_aln.tStart;
                    int new_qStartPos = qStartPos;
                    int new_tStartPos = tStartPos;
                    
                    if (qStartPos == -1 || tStartPos == -1 || alnLen < 3) {
                        continue;
                    }

                    // Find seed point (3 consecutive matches)
                    bool foundConsecutiveMatchSeed = false;
                    for (int blockIdx = 0; blockIdx <= alnLen - 3; ++blockIdx) {
                        int qpos = qStartPos + blockIdx;
                        int dbpos = tStartPos + blockIdx;
                        
                        if (querySequence[qpos] == targetSequence[dbpos] &&
                            querySequence[qpos + 1] == targetSequence[dbpos + 1] &&
                            querySequence[qpos + 2] == targetSequence[dbpos + 2]) {
                            new_qStartPos = qpos + 1; 
                            new_tStartPos = dbpos + 1;
                            foundConsecutiveMatchSeed = true;
                            break;
                        }
                    }
                    
                    if (foundConsecutiveMatchSeed) {
                        std::string backtrace;
                        // Option 1: block alignment (forward + backward)
                        // s_align alignment = blockaligner.align(&target, new_qStartPos, new_tStartPos, backtrace, x_drop, par.covThr, par.covMode);
                        // Option 2: banded alignment (bidirectional)
                        s_align alignment = blockaligner.bandedalign(&target, new_qStartPos, new_tStartPos, backtrace, x_drop, par.covThr, par.covMode);
                        unsigned int alnLength = backtrace.size();
                        double seqId = Util::computeSeqId(par.seqIdMode, alignment.identicalAACnt, query.L, targetLen, alnLength);
                        Matcher::result_t res = Matcher::result_t(targetKey, alignment.score1, alignment.qCov, alignment.tCov, seqId, alignment.evalue, alnLength,
                                                    alignment.qStartPos1, alignment.qEndPos1, query.L, alignment.dbStartPos1, alignment.dbEndPos1, targetLen, backtrace);
                        if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            resultBuffer[bufferIdx].memberIds.push_back(targetId);
                        }
                    }
                }
            } // End of parallel for loop
            
            // Wait for all threads to complete their work
#pragma omp barrier
            
            // Master thread applies results sequentially to maintain priority
#pragma omp master
            {   
                std::vector<unsigned int> keys;
                for (size_t c = 0; c < availChunkSize; ++c) {
                    keys.clear();
                    // Skip invalid or already processed items
                    if (!resultBuffer[c].isValid) {
                        continue;
                    }
                    unsigned int representative = resultBuffer[c].memberIds[0];
                    if (resultBuffer[c].memberIds.size() <= 1 || assignedcluster[representative] != UINT_MAX) {
                        continue;
                    }
                    // Assign all members to this cluster
                    for (size_t id = 0; id < resultBuffer[c].memberIds.size(); id++) {
                        unsigned int memberId = resultBuffer[c].memberIds[id];
                        if (assignedcluster[memberId] == UINT_MAX) {
                            keys.push_back(memberId);
                        }
                    }
                    if (keys.size() <= 1) {
                        continue;
                    }
                    // Assign cluster
                    for (size_t id = 0; id < keys.size(); id++) {
                        assignedcluster[keys[id]] = representative;
                    }
                }
            }
            
            // Wait for master to complete before moving to next chunk
#pragma omp barrier
        } // End of chunk loop
    } // End of parallel region
    
    // Correct edges that are not assigned properly
    for (size_t id = 0; id < dbSize; ++id) {
        if (assignedcluster[id] == UINT_MAX) {
            assignedcluster[id] = id;
        }
    }
    
   
    std::pair<unsigned int, unsigned int> * assignment = new std::pair<unsigned int, unsigned int> [dbSize];
#pragma omp parallel
    {

#pragma omp for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            if (assignedcluster[i] == UINT_MAX) {
                Debug(Debug::ERROR) << "there must be an error: "<< "\t" << i <<
                                    "\tis not assigned to a cluster\n";
                continue;
            }

            assignment[i].first = seqDbr->getDbKey(assignedcluster[i]);
            assignment[i].second = seqDbr->getDbKey(i);
        }
    }
    SORT_PARALLEL(assignment,assignment+dbSize);

    size_t cluNum = (dbSize > 0) ? 1 : 0; // gg // why?
    for(size_t i = 1; i < dbSize; i++){
        cluNum += (assignment[i].first != assignment[i-1].first);
    }

    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << cluNum << "\n\n";
    // write Data
    writeData(&resultWriter, assignment, dbSize);
    // Clean up (add your result writing code here)
    delete[] assignedcluster;
    delete subMat;
    seqDbr->close();
    delete seqDbr;
    
    return 0;
}

int alignblock(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;
    timer.reset();
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    int dbtype = alnDbr.getDbtype();

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, dbtype);
    resultWriter.open();

    int status = doAlign2clust(par, resultWriter, alnDbr, 0, alnDbr.getSize());
    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";
    resultWriter.close();
    alnDbr.close();
    return status;
}