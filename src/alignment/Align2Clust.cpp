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


struct ClusterEntry {
    unsigned int qKey; // rep key
    unsigned int prefSize;
    size_t clusterEntryPosition;
};

struct ClusterSize {
    unsigned int qKey;
    unsigned int cluSize;

    static bool compareByClusterSizeNqKey(const ClusterSize& a, const ClusterSize& b) {
        if (a.cluSize > b.cluSize) {
            return true;
        }
        if (a.cluSize < b.cluSize) {
            return false;
        }
        return a.qKey < b.qKey;
    }
};

struct Element {
    unsigned int key;
    unsigned short diagonal;
};

struct AssignClusterThread {
    size_t priority;      // elementLookupTable index (priority by cluster size)
    unsigned int qKey;
    std::vector<unsigned int> keys;
    AssignClusterThread(): priority(SIZE_MAX),qKey(UINT_MAX),keys() {}
    AssignClusterThread(size_t priority, unsigned int qKey, size_t elementSize)
        : priority(priority), qKey(qKey), keys(elementSize, UINT_MAX) {}
};

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
    // std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
    //                                 ? getCovSeqidQscPercMinDiag()
    //                                 : getCovSeqidQscPercMinDiagTargetCov();
                                    
    // scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    std::cout << "Score per column threshold for filtering: " << scorePerColThr << "\n";
    EvalueComputation evaluer(seqDbr->getAminoAcidDBSize(), subMat);

    // setup linear data
    unsigned int *assignedcluster = new(std::nothrow) unsigned int[dbSize];
    Util::checkAllocation(assignedcluster, "Can not allocate assignedcluster memory in Align2Clust");
    std::fill_n(assignedcluster, dbSize, UINT_MAX);

    // setup mutlithreading
    std::vector<BlockAligner *> blockAligner;
    std::vector<Sequence *> targetSeq;
    std::vector<Matcher *> matcher;
    blockAligner.resize(par.threads);
    targetSeq.resize(par.threads);
    matcher.resize(par.threads);

    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
        blockAligner[thread_idx] = new BlockAligner(par.maxSeqLen, MIN_SIZE, MAX_SIZE, -gapOpen, -gapExtend, par.compBiasCorrection, par.compBiasCorrectionScale, par, Parameters::DBTYPE_AMINO_ACIDS);
        targetSeq[thread_idx] = new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        matcher[thread_idx] = new Matcher(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
#endif
    }


    // setup query
    Sequence querySeq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);

    const size_t flushSize = 1000000;
    Debug::Progress progress(dbSize);
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbSize)/static_cast<double>(flushSize)));

    std::vector<std::pair<unsigned int, unsigned short>> dbKeysNdiagonal;
    std::vector<unsigned int> assignedMemKeys;
    for (size_t i = 0; i < alnDbr.getSize(); i++) {
        const unsigned int queryKey = seqDbr->getDbKey(i);
        const size_t alnId = alnDbr.getId(queryKey);
        char *data = alnDbr.getData(alnId, 0);
        // only load query data if data != \0
        size_t representative;
        // if (*data != '\0') {
            size_t queryId = seqDbr->getId(queryKey);
            if (assignedcluster[queryId] != UINT_MAX) {
                continue;
            }
            representative = queryId;
            char*  querySequence = seqDbr->getData(queryId, 0);
            size_t queryLen = seqDbr->getSeqLen(queryId);
            querySeq.mapSequence(queryId, queryKey, querySequence, queryLen);
        // }

        dbKeysNdiagonal.clear();
        assignedMemKeys.clear();
        while (*data != '\0') {
            // char dbKeyBuffer[256];
            // Util::parseKey(data, dbKeyBuffer);
            // const unsigned int dbKey = static_cast<unsigned int>(strtoul(dbKeyBuffer, NULL, 10));
            // dbKeys.push_back(dbKey);
            hit_t result = QueryMatcher::parsePrefilterHit(data);
            if (assignedcluster[seqDbr->getId(result.seqId)] != UINT_MAX) {
                data = Util::skipLine(data);
                continue;
            }
            dbKeysNdiagonal.push_back(std::make_pair(result.seqId, result.diagonal));
            data = Util::skipLine(data);
        }
        if (dbKeysNdiagonal.size() == 1){
            //singleton
            assignedcluster[representative] = representative;
            continue;
        }

#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
            blockAligner[thread_idx]->initQuery(querySequence, querySeq.numSequence, querySeq.L, Parameters::DBTYPE_AMINO_ACIDS);
            matcher[thread_idx]->initQuery(&querySeq);
#endif
        }

#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            #pragma omp for schedule(dynamic, 1)
            for (size_t id = 0; id < dbKeysNdiagonal.size(); id++) {
                const unsigned int targetKey = dbKeysNdiagonal[id].first;
                const unsigned short diagonal = dbKeysNdiagonal[id].second;
                const size_t targetId = seqDbr->getId(targetKey);
                // if query and target are identical
                const bool isIdentity = (queryKey == targetKey && par.includeIdentity) ? true : false;
                if (queryKey == targetKey) {
                    assignedcluster[targetId] = representative;
                    continue;
                }
                if (isIdentity) {
                    EXIT(EXIT_FAILURE);
                }
                char* targetSequence = seqDbr->getData(targetId, thread_idx);
                size_t targetLen = seqDbr->getSeqLen(targetId);
                targetSeq[thread_idx]->mapSequence(targetId, targetKey, targetSequence, targetLen);

                bool canBeClustered = false;

                // 1. run ungapped alignment
                BlockAligner::LocalAln ungapped_aln = blockAligner[thread_idx]->ungappedAlign(targetSequence, targetLen, diagonal, fastMatrix.matrix); 
                // 1.2. parse ungapped alignment
                int distance = ungapped_aln.score;
                int qStartPos = ungapped_aln.qStart;
                int qEndPos = ungapped_aln.qEnd;
                int tStartPos = ungapped_aln.tStart;
                int tEndPos = ungapped_aln.tEnd;

                //check ungapped criteria
                double evalue = evaluer.computeEvalue(distance, querySeq.L);
                int bitScore = static_cast<int>(evaluer.computeBitScore(distance) + 0.5);
                int alnLen = (qEndPos - qStartPos) + 1;
                // debug
                float queryCov = SmithWaterman::computeCov(qStartPos, qEndPos , querySeq.L);
                float targetCov = SmithWaterman::computeCov(tStartPos, tEndPos , targetLen);

                float currScorePerCol = static_cast<float>(distance) / static_cast<float>(ungapped_aln.diagonalLen);
                bool hasEvalue = (evalue <= par.evalThr);
                bool hasAlnLen = (alnLen >= par.alnLenThr);
                bool hasCov = Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov);
                double seqId = 0;

                const char* querySequence = querySeq.getSeqData();
                if (hasEvalue) {
                    int idCnt = 0;
                    for (int q = qStartPos; q <= qEndPos; q++) {
                        char qLetter = querySequence[q] & static_cast<unsigned char>(~0x20);
                        char tLetter = targetSequence[tStartPos + (q - qStartPos)] & static_cast<unsigned char>(~0x20);
                        idCnt += (qLetter == tLetter) ? 1 : 0;
                    }
                    seqId = Util::computeSeqId(par.seqIdMode, idCnt, querySeq.L, targetLen, alnLen);
                }
                
                bool hasToFilter = (par.filterHits == true && currScorePerCol >= scorePerColThr);
                hasToFilter = 0;
                bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                
                
                if (hasToFilter || (hasAlnLen && hasCov && hasSeqId && hasEvalue)) {
                    canBeClustered = true;
                } else {
                    // Run gapped alignment
                    std::string backtrace;
                    int32_t x_drop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());
                    int new_qStartPos = qStartPos;
                    int new_tStartPos = tStartPos;
                    
                    bool foundMatch = false;
                    unsigned char *targetNumSequence = targetSeq[thread_idx]->numSequence;
            
                    for (int j=0; j < alnLen; ++j){
                        int qpos = qStartPos + j;
                        int dbpos = tStartPos + j;
                        if (querySeq.numSequence[qpos] == targetNumSequence[dbpos]) {
                            new_qStartPos = qpos;
                            new_tStartPos = dbpos;
                            foundMatch = true;
                            break;
                        }
                    }
                    
                    if(foundMatch) {
                        //temporary
                        s_align alignment = blockAligner[thread_idx]->align(targetSequence, targetNumSequence, targetLen, new_qStartPos, new_tStartPos, backtrace, &evaluer, x_drop);
                        unsigned int alnLength = backtrace.size(); // Is it correct?
                        double seqId = Util::computeSeqId(par.seqIdMode, alignment.identicalAACnt, querySeq.L, targetLen, alnLength);
                        Matcher::result_t res = Matcher::result_t(targetKey, alignment.score1, alignment.qCov, alignment.tCov, seqId, alignment.evalue, alnLen,
                                                    alignment.qStartPos1, alignment.qEndPos1, querySeq.L, alignment.dbStartPos1, alignment.dbEndPos1, targetLen, backtrace);
                        if (Alignment::checkCriteria(res, 0, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            canBeClustered = true;
                        }
                    } 
                }

                if (canBeClustered) {
                    dbKeysNdiagonal[id].second = USHRT_MAX; // mark as clustered
                    assignedcluster[targetId] = representative; // assign representative
                }
            } 
        }
        // Expand cluster
        std::set<unsigned int> newlyAssignedMemIds;
        #pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            std::vector<unsigned int> localAssignedMemIds;

            #pragma omp for schedule(dynamic, 1)
            for (size_t id = 0; id < dbKeysNdiagonal.size(); id++) {
                if (dbKeysNdiagonal[id].second != USHRT_MAX) continue; // not clustered
                const unsigned int assignedMemKey = dbKeysNdiagonal[id].first;
                const size_t alnId = alnDbr.getId(assignedMemKey);
                char* dataExpand = alnDbr.getData(alnId, thread_idx);

                while (*dataExpand != '\0') {
                    hit_t result = QueryMatcher::parsePrefilterHit(dataExpand);
                    
                    // if it's already assigned, skip
                    if (assignedcluster[seqDbr->getId(result.seqId)] != UINT_MAX) {
                        dataExpand = Util::skipLine(dataExpand);
                        continue;
                    }
                    const unsigned int targetKey = result.seqId;
                    const size_t targetId = seqDbr->getId(targetKey);

                    // if query and target are identical
                    const bool isIdentity = (queryKey == targetKey && par.includeIdentity) ? true : false;
                    if (isIdentity) {
                        continue;
                    }
                    char* targetSequence = seqDbr->getData(targetId, thread_idx);
                    size_t targetLen = seqDbr->getSeqLen(targetId);
                    targetSeq[thread_idx]->mapSequence(targetId, targetKey, targetSequence, targetLen);

                    // run Alignment
                    Matcher::result_t resultSW = matcher[thread_idx]->getSWResult(targetSeq[thread_idx], INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, 0);
                    if (Alignment::checkCriteria(resultSW, 0, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        localAssignedMemIds.push_back(targetId);
                    }
                    dataExpand = Util::skipLine(dataExpand);
                }

            }

            // merge local assigned mem ids to global set
            #pragma omp critical
            {
                for (auto& id : localAssignedMemIds) {
                    newlyAssignedMemIds.insert(id);
                }
            }

            // assign representative to newly assigned members
            for (auto& id : newlyAssignedMemIds) {
                if(assignedcluster[id] == UINT_MAX){
                    assignedcluster[id] = representative;
                }
            }
        }
    }
    
    // correct edges that are not assigned properly
    for (size_t id = 0; id < dbSize; ++id) {
        // check if the assigned clusterid is a rep. sequence
        // if not, make it a rep. seq. again
        if(assignedcluster[id] == UINT_MAX){
            std::cout << "unassigned id found: " << id << "\n";
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

            // assignment[i].first = seseqDbr->getDbKey(assignedcluster[i]);
            // assignment[i].second = seseqDbr->getDbKey(i);
            assignment[i].first = assignedcluster[i];
            assignment[i].second = i;
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

    // delete
    delete[] assignedcluster;
    delete[] assignment;

    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
    delete subMat;
    return 0;
}

int align2clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    int dbtype = alnDbr.getDbtype();

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, dbtype);
    resultWriter.open();

    int status = doAlign2clust(par, resultWriter, alnDbr, 0, alnDbr.getSize());
    resultWriter.close();
    alnDbr.close();
    return status;
}