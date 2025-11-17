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
#include "EvalueComputation.h"
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
    tdbr = tDbrIdx->sequenceReader;
    bool sameQTDB = (par.db2.compare(par.db1) == 0);
    if (sameQTDB == true) {
        qDbrIdx = tDbrIdx;
        qdbr = tdbr;
    } else {
        // open the sequence, prefiltering and output databases
        qDbrIdx = new IndexReader(par.db1, par.threads,  IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        qdbr = qDbrIdx->sequenceReader;
    }

    if(resultDbr.isSortedByOffset() && qdbr->isSortedByOffset()){
        qdbr->setSequentialAdvice();
    }

    BaseMatrix *subMat;
    subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    float scorePerColThr = 0.0;
    std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
                                    
    scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    Debug(Debug::INFO) << "Score per column threshold for filtering: " << scorePerColThr << "\n";
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
            // Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
            char buffer[1024 + 32768*4];
            std::string resultBuffer;
            resultBuffer.reserve(1000000);
            Sequence query(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
            BlockAligner blockaligner(Parameters::DBTYPE_AMINO_ACIDS, par.maxSeqLen, subMat, &fastMatrix, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
            
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
                    const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameQTDB)) ? true : false;
                    const char* targetSeq = tdbr->getData(targetId, thread_idx);
                    size_t targetLen = tdbr->getSeqLen(targetId);
                    target.mapSequence(targetId, targetKey, targetSeq, targetLen);
                    
                    double seqId = 0.0;
                    bool hasAlnLen = false;
                    bool hasCov = false;
                    bool hasSeqId = false;
                    std::string backtrace = "";
                    // 0. run hamming alignment
                    if (par.skipHamming == false){
                        BlockAligner::UngappedAln_res hamming_aln = blockaligner.hammingDistance(&target, prefRes[element].diagonal); 
                        double seqId = Util::computeSeqId(par.seqIdMode, static_cast<float>(hamming_aln.score), query.L, targetLen, hamming_aln.diagonalLen);
                        bool hasAlnLen = (hamming_aln.alnLen >= par.alnLenThr);
                        float hammingCovThr = std::max(0.5f, par.covThr);
                        float hammingSeqIdThr = std::max(0.5f, par.seqIdThr);
                        bool hasCov = Util::hasCoverage(hammingCovThr, par.covMode, hamming_aln.qcov, hamming_aln.tcov);
                        bool hasSeqId = seqId >= (hammingSeqIdThr - std::numeric_limits<float>::epsilon());
                        if (hasAlnLen && hasCov && hasSeqId) {
                            Matcher::result_t result;
                            result = Matcher::result_t(targetKey, hamming_aln.bitScore, hamming_aln.qcov, hamming_aln.tcov, seqId, hamming_aln.eval, hamming_aln.alnLen,
                                                        hamming_aln.qStart, hamming_aln.qEnd, query.L, hamming_aln.tStart, hamming_aln.tEnd, target.L, backtrace);
                            alnResults.emplace_back(result);
                            continue;
                        }
                    }
                    
                    
                    
                    // 1. run ungapped alignment
                    BlockAligner::UngappedAln_res ungapped_aln = blockaligner.ungappedAlign(&target, prefRes[element].diagonal); 
                    //check ungapped criteria
                    bool hasEvalue = (ungapped_aln.eval <= par.evalThr);
                    hasAlnLen = (ungapped_aln.alnLen >= par.alnLenThr);
                    hasCov = Util::hasCoverage(par.covThr, par.covMode, ungapped_aln.qcov, ungapped_aln.tcov);
                    seqId = 0;
                    if (hasEvalue) {    
                        int idCnt = 0;
                        for (int q = ungapped_aln.qStart; q <= ungapped_aln.qEnd; q++) {
                            char qLetter = querySeq[q] & static_cast<unsigned char>(~0x20);
                            char tLetter = targetSeq[ungapped_aln.tStart + (q - ungapped_aln.qStart)] & static_cast<unsigned char>(~0x20);
                            idCnt += (qLetter == tLetter) ? 1 : 0;
                        }
                        seqId = Util::computeSeqId(par.seqIdMode, idCnt, query.L, target.L, ungapped_aln.alnLen);
                    }
                    char *end = Itoa::i32toa_sse2(ungapped_aln.alnLen, buffer);
                    size_t len = end - buffer;
                    backtrace = "";
                    if (par.addBacktrace) {
                        backtrace=std::string(buffer, len - 1);
                        backtrace.push_back('M');
                    }
                    if (isIdentity) {
                        // set coverage and seqid of identity
                        ungapped_aln.qcov = 1.0f;
                        ungapped_aln.tcov = 1.0f;
                        seqId = 1.0f;
                    }
                    hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    
                    if (isIdentity || (hasAlnLen && hasCov && hasSeqId && hasEvalue)) {
                        Matcher::result_t result;
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

                    int qStartPos = ungapped_aln.qStart;
                    int tStartPos = ungapped_aln.tStart;

                    int new_qStartPos = qStartPos;
                    int new_tStartPos = tStartPos;
                    if (qStartPos == -1 || tStartPos == -1 || alnLen < 3) {
                        continue;
                    }

                    bool foundConsecutiveMatchSeed = false;
    
                    for (int blockIdx = 0; blockIdx <= alnLen - 3; ++blockIdx) {
                        int qpos = qStartPos + blockIdx;
                        int dbpos = tStartPos+ blockIdx;
    
    
                        if (querySeq[qpos] == targetSeq[dbpos] &&
                            querySeq[qpos + 1] == targetSeq[dbpos + 1] &&
                            querySeq[qpos + 2] == targetSeq[dbpos + 2]) {
                            
                            // Found 3 consecutive matches
                            new_qStartPos = qpos + 1; 
                            new_tStartPos = dbpos + 1;
                            foundConsecutiveMatchSeed = true;
                            break;
                        }
                    }
            
                    if (foundConsecutiveMatchSeed){
                        std::string backtrace;
                        // Option1. block alignment (forward + backward)
                        // s_align alignment = blockaligner.align(&target, new_qStartPos, new_tStartPos, backtrace, x_drop, par.covThr, par.covMode);
                        // Option2. banded alignment (bidirectional)
                        s_align alignment = blockaligner.bandedalign(&target, new_qStartPos, new_tStartPos, backtrace, x_drop, par.covThr, par.covMode);
                        unsigned int alnLength = backtrace.size();
                        double seqId = Util::computeSeqId(par.seqIdMode, alignment.identicalAACnt, query.L, targetLen, alnLength);
                        Matcher::result_t res = Matcher::result_t(targetKey, alignment.score1, alignment.qCov, alignment.tCov, seqId, alignment.evalue, alnLength,
                                                    alignment.qStartPos1, alignment.qEndPos1, query.L, alignment.dbStartPos1, alignment.dbEndPos1, targetLen, backtrace);
                        if (Alignment::checkCriteria(res,isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            alnResults.emplace_back(res);
                        }
                    } 
                }

                if (par.sortResults > 0 && alnResults.size() > 1) {
                    SORT_SERIAL(alnResults.begin(), alnResults.end(), Matcher::compareHits);
                }

                for (size_t i = 0; i < alnResults.size(); ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, alnResults[i], par.addBacktrace); 
                    resultBuffer.append(buffer, len);
                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
                resultBuffer.clear();
                alnResults.clear();
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
    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
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