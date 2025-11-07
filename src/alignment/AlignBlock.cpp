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


#define MAX_SIZE 256
#define MIN_SIZE 32

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

int doalignblock(Parameters &par,
                DBWriter &resultWriter,
                DBReader<unsigned int> &resultDbr,
            const size_t dbFrom, const size_t dbSize) {
    IndexReader * qDbrIdx = NULL;
    DBReader<unsigned int> * qdbr = NULL;
    DBReader<unsigned int> * tdbr = NULL;
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader * tDbrIdx = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES,   (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    int querySeqType = 0;
    tdbr = tDbrIdx->sequenceReader;
    int targetSeqType = tDbrIdx->getDbtype();
    bool sameQTDB = (par.db2.compare(par.db1) == 0);
    if (sameQTDB == true) {
        qDbrIdx = tDbrIdx;
        qdbr = tdbr;
        querySeqType = targetSeqType;
    } else {
        // open the sequence, prefiltering and output databases
        qDbrIdx = new IndexReader(par.db1, par.threads,  IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        qdbr = qDbrIdx->sequenceReader;
        querySeqType = qdbr->getDbtype();
    }

    if(resultDbr.isSortedByOffset() && qdbr->isSortedByOffset()){
        qdbr->setSequentialAdvice();
    }

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
    std::cout << "Score per column threshold for filtering: " << scorePerColThr << "\n";
    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), subMat);

    int32_t x_drop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());

    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 100000000;
    if (totalMemory > resultDbr.getTotalDataSize()) {
        flushSize = resultDbr.getSize();
    }
    
    size_t iterations = 1;
    if(flushSize > 0){
        iterations = static_cast<int>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));
    }

    for (size_t iter = 0; iter < iterations; iter++) {
        size_t start = dbFrom + (iter * flushSize);
        size_t bucketSize = std::min(dbSize - (iter * flushSize), flushSize);
        Debug::Progress progress(bucketSize);

#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char buffer[1024 + 32768*4];
            std::string resultBuffer;
            // Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
            Sequence query(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
            BlockAligner blockaligner(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
            
            std::vector<Matcher::result_t> alnResults;
            alnResults.reserve(300);
#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                
                char *data = resultDbr.getData(id, thread_idx);
                size_t queryKey = resultDbr.getDbKey(id);
                const char* querySeq = NULL;
                if(*data !=  '\0'){
                    const size_t queryId = qdbr->getId(queryKey);
                    querySeq = qdbr->getData(queryId, thread_idx);
                    size_t queryLen = static_cast<int>(qdbr->getSeqLen(queryId));
                    query.mapSequence(queryId, queryKey, querySeq, queryLen);
                    blockaligner.initQuery(&query);
                }

                std::vector<hit_t> prefRes = QueryMatcher::parsePrefilterHits(data);
                for (size_t element = 0; element < prefRes.size(); element++) {
                    const unsigned int targetKey = tdbr->getDbKey(prefRes[element].seqId);
                    const unsigned int targetId = tdbr->getId(targetKey);
                    const char* targetSeq = tdbr->getData(targetId, thread_idx);
                    size_t targetLen = tdbr->getSeqLen(targetId);
                    target.mapSequence(targetId, targetKey, targetSeq, targetLen);

                    // 1. run ungapped alignment
                    BlockAligner::UngappedAln_res ungapped_aln = blockaligner.ungappedAlign(&target, prefRes[element].diagonal); 
                    //check ungapped criteria
                    bool hasEvalue = (ungapped_aln.eval <= par.evalThr);
                    bool hasAlnLen = (ungapped_aln.alnLen >= par.alnLenThr);
                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, ungapped_aln.qcov, ungapped_aln.tcov);
                    double seqId = 0;
                    if (hasEvalue) {    
                        int idCnt = 0;
                        for (int q = ungapped_aln.qStart; q <= ungapped_aln.qEnd; q++) {
                            char qLetter = querySeq[q] & static_cast<unsigned char>(~0x20);
                            char tLetter = targetSeq[ungapped_aln.tStart + (q - ungapped_aln.qStart)] & static_cast<unsigned char>(~0x20);
                            idCnt += (qLetter == tLetter) ? 1 : 0;
                        }
                        seqId = Util::computeSeqId(par.seqIdMode, idCnt, query.L, target.L, ungapped_aln.alnLen);
                    }
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    
                    Matcher::result_t result;
                    std::string backtrace = "";
                    if (hasAlnLen && hasCov && hasSeqId && hasEvalue) {
                        result = Matcher::result_t(targetKey, ungapped_aln.bitScore, ungapped_aln.qcov, ungapped_aln.tcov, seqId, ungapped_aln.eval, ungapped_aln.alnLen,
                                                       ungapped_aln.qStart, ungapped_aln.qEnd, query.L, ungapped_aln.tStart, ungapped_aln.tEnd, target.L, backtrace);
                        alnResults.emplace_back(result);
                        continue;
                    }
                        
                    // If ungapped criteria not met
                    // 2. run gapped alignment
                    // 2-1. Check if we need to run gapped alignment
                    float currScorePerCol = static_cast<float>(ungapped_aln.score) / static_cast<float>(ungapped_aln.diagonalLen);
                    if (currScorePerCol < scorePerColThr) { // if fail score per column threshold, skip
                        continue;
                    } 
                    // 2-2. Run gapped alignment
                    // find seed point from mid
                    int alnLen = ungapped_aln.alnLen;
                    int mid = alnLen / 2;

                    int qStartPos = ungapped_aln.qStart;
                    int tStartPos = ungapped_aln.tStart;

                    int new_qStartPos = qStartPos;
                    int new_tStartPos = tStartPos;

                    bool foundMatch = false;

                    for (int offset = 0; offset < alnLen; ++offset) {
                        // mid - offset (left)
                        int j_left = mid - offset;
                        if (j_left >= 0) {
                            int qpos = qStartPos + j_left;
                            int dbpos = tStartPos + j_left;
                            if (querySeq[qpos] == targetSeq[dbpos]) {
                                new_qStartPos = qpos;
                                new_tStartPos = dbpos;
                                foundMatch = true;
                                break;
                            }
                        }

                        // mid + offset (right)
                        int j_right = mid + offset;
                        if (j_right < alnLen) {
                            int qpos = qStartPos + j_right;
                            int dbpos = tStartPos + j_right;
                            if (querySeq[qpos] == targetSeq[dbpos]) {
                                new_qStartPos = qpos;
                                new_tStartPos = dbpos;
                                foundMatch = true;
                                break;
                            }
                        }
                    }
                
                    if(foundMatch) {
                        s_align alignment = blockaligner.align(&target, new_qStartPos, new_tStartPos, backtrace, x_drop);
                        unsigned int alnLength = backtrace.size(); // Is it correct?
                        double seqId = Util::computeSeqId(par.seqIdMode, alignment.identicalAACnt, query.L, targetLen, alnLength);
                        Matcher::result_t res = Matcher::result_t(targetKey, alignment.score1, alignment.qCov, alignment.tCov, seqId, alignment.evalue, alnLen,
                                                    alignment.qStartPos1, alignment.qEndPos1, query.L, alignment.dbStartPos1, alignment.dbEndPos1, targetLen, backtrace);
                        if (Alignment::checkCriteria(res, 0, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            alnResults.emplace_back(res);
                        }
                    } 
                }

                if (par.sortResults > 0 && alnResults.size() > 1) {
                    SORT_SERIAL(alnResults.begin(), alnResults.end(), Matcher::compareHits);
                }
                for (size_t i = 0; i < alnResults.size(); ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, alnResults[i], false, false); // add backtrace false   // gyuri
                    resultBuffer.append(buffer, len);
                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
                resultBuffer.clear();
            }
        }
    }
     if (tDbrIdx != NULL) {
        delete tDbrIdx;
    }

    if (sameQTDB == false) {
        if(qDbrIdx != NULL){
            delete qDbrIdx;
        }
    }
    // delete
    delete subMat;
    return 0;
}

int alignblock(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;
    timer.reset();
    DBReader<unsigned int> resultDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    int dbtype = resultDbr.getDbtype();

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, dbtype);
    resultWriter.open();

    int status = doalignblock(par, resultWriter, resultDbr, 0, resultDbr.getSize());
    Debug(Debug::INFO) << "Time for run alignblock: " << timer.lap() << " sec\n";
    resultWriter.close();
    resultDbr.close();
    return status;
}