#include "BlockAligner.h"
#include "Sequence.h"
#include "Util.h"
#include "Parameters.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"

#define MAX_SIZE 4096
#define MIN_SIZE 32

BlockAligner::BlockAligner(
    int dbtype,
    size_t maxSequenceLength,
    BaseMatrix *m, SubstitutionMatrix::FastMatrix* fastMatrix, EvalueComputation * evaluer,
    bool compBiasCorrection, float compBiasCorrectionScale,
    int8_t gapOpen, int8_t gapExtend
) : 
    gaps({gapOpen, gapExtend}),
    compBiasCorrection(compBiasCorrection),
    compBiasCorrectionScale(compBiasCorrectionScale),
    dbtype(dbtype),
    subMat((SubstitutionMatrix*) m),
    fastMatrix(fastMatrix),
    evaluer(evaluer) {
    range={MIN_SIZE, MAX_SIZE};
    query = block_new_padded_aa(maxSequenceLength, range.max);
    queryBias = block_new_pos_bias(maxSequenceLength, range.max);
    
    
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) == false) {
        querySeqType = Parameters::DBTYPE_AMINO_ACIDS;
        target = block_new_padded_aa(maxSequenceLength, range.max);
        targetBias = block_new_pos_bias(maxSequenceLength, range.max);
        matrix = block_new_simple_aamatrix(1, -1);
        for (int aa1 = 0; aa1 < subMat->alphabetSize; aa1++) {
			for (int aa2 = 0; aa2 < subMat->alphabetSize; aa2++) {
				// instead of num2aa, use aa directly
				block_set_aamatrix_num(matrix, aa1, aa2,
								subMat->subMatrix[aa1][aa2]);
			}
		}
    } else {
        querySeqType = Parameters::DBTYPE_HMM_PROFILE;
        bProfile = block_new_aaprofile(maxSequenceLength, range.max, gaps.extend);
    }
    blockTrace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, range.max);
    blockNoTrace = block_new_aa_xdrop(maxSequenceLength, maxSequenceLength, range.max);
    cigar = block_new_cigar(maxSequenceLength, maxSequenceLength);
    sAlnCigar = new uint32_t[maxSequenceLength];
    queryCompBias = new int16_t[maxSequenceLength];
    targetCompBias = new int16_t[maxSequenceLength];
    //set targetCompBias to 0
    memset(targetCompBias, 0, maxSequenceLength * sizeof(int16_t));
    tmpCompBias   = new float[maxSequenceLength];
}

BlockAligner::~BlockAligner() {
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(blockTrace);
    block_free_aa_xdrop(blockNoTrace);
    block_free_padded_aa(query);
    block_free_pos_bias(queryBias);
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) == false) {
        block_free_padded_aa(target);
        block_free_pos_bias(targetBias);
        block_free_aamatrix(matrix);
    } else {
        block_free_aaprofile(bProfile);
    }
    delete[] sAlnCigar;
    delete[] queryCompBias;
    delete[] targetCompBias;
    delete[] tmpCompBias;

}

void BlockAligner::initQuery(Sequence* query){
    currentQuery = query;
    this->querySeq = query->getSeqData();
    this->queryNumSeq = query->numSequence;
    this->queryLength = query->L;

    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)&& compBiasCorrection){
        SubstitutionMatrix::calcLocalAaBiasCorrection(subMat, this->queryNumSeq, this->queryLength, tmpCompBias, compBiasCorrectionScale);
        for (int i =0; i < this->queryLength; i++) { 
			// queryCompBias[i] = (int8_t) (tmpCompBias[i] < 0.0)? tmpCompBias[i] - 0.5: tmpCompBias[i] + 0.5; //why .5 offset?
            queryCompBias[i] = (int8_t) (tmpCompBias[i]);
		}
    }
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlign(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln;

    AlignResult res;
    res.score = -1000000000;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;

    // forwards alignment starting at (qIdx, tIdx)
    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(qNum + qIdx), qLen - qIdx, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(tNum + tIdx), tLen - tIdx, range.max);
    
    // PosBias
    block_set_pos_bias(queryBias, queryCompBias + qIdx, qLen - qIdx);
    block_set_pos_bias(targetBias, targetCompBias + tIdx, tLen - tIdx);

    // t_bias = block_new_pos_bias(tLen - tIdx, range.max); // or comment it out since it is full of 0
    block_align_aa_xdrop_posbias(blockNoTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); // forward with no trace

    res = block_res_aa_xdrop(blockNoTrace);

    int qEnd = qIdx + res.query_idx - 1;
    int tEnd = tIdx + res.reference_idx - 1;
    // float maxQueryCov = SmithWaterman::computeCov(0, qEnd , qLen);
    // float maxTargetCov = SmithWaterman::computeCov(0, res_aln.tEnd , tLen);
    // bool hasCov = Util::hasCoverage(par.covThr, par.covMode, maxQueryCov, maxTargetCov);
    // if (!hasCov) {
    //     local_aln.evalue = -1.0f; // this should avoid that the hit is added
    //     local_aln.score1 = -1.0f;
    //     local_aln.qCov = 0.0f;
    //     local_aln.tCov = 0.0f;
    //     return local_aln;
    // }
    
    // Debug(Debug::ERROR) << "end position: " << res_aln.qEnd << "\t" << res_aln.tEnd << "\n";
    // if (res_aln.qEnd <=0 || res_aln.tEnd <=0) {
    if (qEnd == SIZE_MAX || tEnd == SIZE_MAX) {
        // Debug(Debug::ERROR) << "wrong end position: " << qEnd << "\t" << tEnd << "\n";
        local_aln.score1 = -1.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }
    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa_numsequence(query, (uint8_t*)qNum, qEnd, range.max);
    block_set_bytes_rev_padded_aa_numsequence(target, (uint8_t*)tNum, tEnd, range.max);
    
    //PosBias
    block_set_rev_pos_bias(queryBias, queryCompBias, qEnd);
    block_set_rev_pos_bias(targetBias, targetCompBias, tEnd);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    block_cigar_eq_aa_trace_xdrop(blockTrace, query, target, res.query_idx, res.reference_idx, cigar);

    int qStart = qEnd - res.query_idx;
    int tStart = tEnd - res.reference_idx;
    int score = res.score;

    float qcov = SmithWaterman::computeCov(qStart, qEnd, qLen);
    float tcov = SmithWaterman::computeCov(tStart, tEnd, tLen);

    double evalue = evaluer->computeEvalue(score, qLen);
    int bitScore = static_cast<int>(evaluer->computeBitScore(score) + 0.5);
    
    local_aln.score1 = bitScore;
    local_aln.qStartPos1 = qStart;
    local_aln.qEndPos1 = qEnd;
    local_aln.dbStartPos1 = tStart;
    local_aln.dbEndPos1 = tEnd;
    local_aln.qCov = qcov;
    local_aln.tCov = tcov;
    local_aln.evalue = evalue;

    return local_aln;
}

BlockAligner::UngappedAln_res align_local_profile(
    BlockHandle blockTrace, BlockHandle blockNoTrace,
    const char* q_str, const size_t qLen, PaddedBytes* a,
    const char* t_str, const size_t tLen, AAProfile* bProfile,
    Gaps gaps, BaseMatrix& subMat,
    int qIdx, int tIdx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    BlockAligner::UngappedAln_res res_aln;
    AlignResult res;

    // forwards alignment starting at (qIdx, tIdx)
    block_set_bytes_padded_aa(a, (uint8_t*)(q_str + qIdx), qLen - qIdx, range.max);

    // assign extra profile columns to 'U', which is unused
    int aa = Sequence::PROFILE_READIN_SIZE;
    uint8_t order[Sequence::PROFILE_READIN_SIZE];
    memset(order, 'U', Sequence::PROFILE_READIN_SIZE);
    memcpy(order, (uint8_t*)subMat.num2aa, Sequence::PROFILE_AA_SIZE);

    block_clear_aaprofile(bProfile, tLen - tIdx, range.max);
    // note: scores are divided by 4 by shifting right by 2
    block_set_all_aaprofile(bProfile, order, aa, (int8_t*)(t_str + tIdx * aa), (tLen - tIdx) * aa, 0, 2);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_xdrop(blockNoTrace, a, bProfile, range, x_drop);
    res = block_res_aa_xdrop(blockNoTrace);

    res_aln.qEnd = qIdx + res.query_idx;
    res_aln.tEnd = tIdx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa(a, (uint8_t*)q_str, res_aln.qEnd, range.max);

    block_clear_aaprofile(bProfile, res_aln.tEnd, range.max);
    block_set_all_rev_aaprofile(bProfile, order, aa, (int8_t*)t_str, res_aln.tEnd * aa, 0, 2);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_trace_xdrop(blockTrace, a, bProfile, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);

    res_aln.qStart = res_aln.qEnd - res.query_idx;
    res_aln.tStart = res_aln.tEnd - res.reference_idx;
    res_aln.score = res.score;
    return res_aln;
}


BlockAligner::UngappedAln_res BlockAligner::ungappedAlign(
    Sequence* target, const unsigned short diagonal)
{
    const char* targetSeq = target->getSeqData();
    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeUngappedAlignment(
                                                        querySeq, queryLength, targetSeq, target->L,
                                                        diagonal, fastMatrix->matrix, Parameters::RESCORE_MODE_ALIGNMENT
                                                    );
    unsigned int distanceToDiagonal = alignment.distToDiagonal;
    int ungappedDiagonal = alignment.diagonal;
    int distance = alignment.score;

    int qUngappedStartPos, qUngappedEndPos, tUngappedStartPos, tUngappedEndPos;
    if (ungappedDiagonal >= 0) {
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        tUngappedStartPos = alignment.startPos;
        tUngappedEndPos = alignment.endPos;
    } else {
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        tUngappedStartPos = alignment.startPos + distanceToDiagonal;
        tUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    //check ungapped criteria
    double evalue = evaluer->computeEvalue(alignment.score, currentQuery->L);
    int bitScore = static_cast<int>(evaluer->computeBitScore(distance) + 0.5);
    int alnLen = (qUngappedEndPos - qUngappedStartPos) + 1;

    float qCov = SmithWaterman::computeCov(qUngappedStartPos, qUngappedEndPos , queryLength);
    float tCov = SmithWaterman::computeCov(tUngappedStartPos, tUngappedEndPos , target->L);
    
    return UngappedAln_res(
        bitScore,
        qCov,
        tCov,
        evalue, 
        alnLen,
        qUngappedStartPos,
        tUngappedStartPos,
        qUngappedEndPos,
        tUngappedEndPos,
        alignment.diagonalLen,
        alignment.score,
        alignment.diagonal
    );
}


s_align
BlockAligner::align(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    //reset sAlnCigar
    memset(sAlnCigar, 0, sizeof(uint32_t) * currentQuery->L);
    s_align local_aln = gappedLocalAlign(currentTarget, qIdx, tIdx, cigar, xdrop);

    int aaIds = 0;
    size_t cigarLength = block_len_cigar(cigar);
    int queryPos = 0;
    int targetPos = 0;
    int queryStartPos = local_aln.qEndPos1;
    int targetStartPos = local_aln.dbEndPos1;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    for (size_t c = 0; c < cigarLength; c++) {
        OpLen o = block_get_cigar(cigar, c);
        char letter = ops_char[o.op];
        switch (letter) {
            case '=':
                aaIds += length;
                // FALLTHROUGH
            case 'X':
                // FALLTHROUGH
            case 'M':
                queryPos += length;
                targetPos += length;
                backtrace.append(length, 'M');
                sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'M');
                break;
            case 'I':
                queryPos += length;
                backtrace.append(length, 'I');
                sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'I');
                break;
            case 'D':
                targetPos += length;
                backtrace.append(length, 'D');
                sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'D');
                break;
        }
        alnLength += length;
        // if(o.op == 1 || o.op == 2){
        //     for(size_t j = 0; j < o.len; j++){
        //         // change traceback with int not char
        //         if(qNum[-queryPos - j + queryStartPos] == tNum[-targetPos - j + targetStartPos]){
        //             aaIds++;
        //         }
        //     }
        //     queryPos += o.len;
        //     targetPos += o.len;
        //     backtrace.append(o.len,'M');
        // }else if(o.op == 4){
        //     queryPos += o.len;
        //     backtrace.append(o.len,'I');
        // }else if(o.op == 5){
        //     targetPos += o.len;
        //     backtrace.append(o.len,'D');
        // }
    }

    std::reverse(backtrace.begin(), backtrace.end());


 
    // // Note: 'M' signals either a match or mismatch
    // char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};

    // int alnLength = 0;
    // size_t cigarLength = block_len_cigar(cigar);
    // size_t aaIds = 0;
    // if (cigarLength > 0) {
    //     int32_t targetPos = 0, queryPos = 0;
    //     for (unsigned long c = 0; c < cigarLength; ++c) {
    //         OpLen o = block_get_cigar(cigar, cigarLength - 1 - c);
    //         char letter = ops_char[o.op];
    //         uint32_t length = o.len;

    //         switch (letter) {
    //             case '=':
    //                 aaIds += length;
    //                 // FALLTHROUGH
    //             case 'X':
    //                 // FALLTHROUGH
    //             case 'M':
    //                 queryPos += length;
    //                 targetPos += length;
    //                 backtrace.append(length, 'M');
    //                 sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'M');
    //                 break;
    //             case 'I':
    //                 queryPos += length;
    //                 backtrace.append(length, 'I');
    //                 sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'I');
    //                 break;
    //             case 'D':
    //                 targetPos += length;
    //                 backtrace.append(length, 'D');
    //                 sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'D');
    //                 break;
    //         }
    //         alnLength += length;
    //     }
    // }

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.cigar = sAlnCigar;
    local_aln.cigarLen = cigarLength;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}



s_align BlockAligner::gappedLocalAlignForward(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln;

    AlignResult res;
    res.score = -1000000000;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;

    // forwards alignment starting at (qIdx, tIdx)
    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(qNum + qIdx), qLen - qIdx, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(tNum + tIdx), tLen - tIdx, range.max);
    
    // PosBias
    block_set_pos_bias(queryBias, queryCompBias + qIdx, qLen - qIdx);
    block_set_pos_bias(targetBias, targetCompBias + tIdx, tLen - tIdx);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); // forward with no trace
    res = block_res_aa_xdrop(blockTrace);
    block_cigar_eq_aa_trace_xdrop(blockTrace, query, target, res.query_idx, res.reference_idx, cigar);

    int qEnd = qIdx + res.query_idx - 1;
    int tEnd = tIdx + res.reference_idx - 1;
    
    // if (res_aln.qEnd <=0 || res_aln.tEnd <=0) {
    if (qEnd == SIZE_MAX || tEnd == SIZE_MAX) {
        Debug(Debug::ERROR) << "wrong end position: " << qEnd << "\t" << tEnd << "\n";
        local_aln.score1 = -1.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }
    
    local_aln.qEndPos1 = qEnd;
    local_aln.dbEndPos1 = tEnd;
    local_aln.score1 = res.score;

    return local_aln;
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlignBackward(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln;

    AlignResult res;
    res.score = -1000000000;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;

    Cigar* cigarBackward = block_new_cigar(res.query_idx, res.reference_idx);
    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa_numsequence(query, (uint8_t*) qNum, qIdx +1, range.max);
    block_set_bytes_rev_padded_aa_numsequence(target, (uint8_t*) tNum, tIdx +1, range.max);
    
    //PosBias
    block_set_rev_pos_bias(queryBias, queryCompBias, qIdx +1 );
    block_set_rev_pos_bias(targetBias, targetCompBias, tIdx +1);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    block_cigar_eq_aa_trace_xdrop(blockTrace, query, target, res.query_idx, res.reference_idx, cigar);

    int qStart = qIdx - res.query_idx +1;
    int tStart = tIdx - res.reference_idx +1;
    int score = res.score;

    local_aln.qStartPos1 = qStart;
    local_aln.dbStartPos1 = tStart;
    local_aln.score1 = res.score;

    return local_aln;
}

s_align
BlockAligner::bandedalignForward(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    memset(sAlnCigar, 0, sizeof(uint32_t) * currentQuery->L);
    //forward
    s_align local_aln = gappedLocalAlignForward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // s_align local_alnBackward = gappedLocalAlignBackward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};
    int alnLength = 0;
    size_t cigarLength = block_len_cigar(cigar);
    size_t aaIds = 0;
    if (cigarLength > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLength; ++c) {
            OpLen o = block_get_cigar(cigar, c);
            char letter = ops_char[o.op];
            uint32_t length = o.len;

            switch (letter) {
                case '=':
                    aaIds += length;
                    // FALLTHROUGH
                case 'X':
                    // FALLTHROUGH
                case 'M':
                    queryPos += length;
                    targetPos += length;
                    backtrace.append(length, 'M');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'M');
                    break;
                case 'I':
                    queryPos += length;
                    backtrace.append(length, 'I');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'I');
                    break;
                case 'D':
                    targetPos += length;
                    backtrace.append(length, 'D');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'D');
                    break;
            }
            alnLength += length;
        }
    }

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.cigar = sAlnCigar;
    local_aln.cigarLen = cigarLength;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}


s_align
BlockAligner::bandedalignBackward(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    //Backward
    memset(sAlnCigar, 0, sizeof(uint32_t) * currentQuery->L);
    s_align local_aln = gappedLocalAlignBackward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // s_align local_alnBackward = gappedLocalAlignBackward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};
    int alnLength = 0;
    size_t cigarLength = block_len_cigar(cigar);
    size_t aaIds = 0;
    if (cigarLength > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLength; ++c) {
            OpLen o = block_get_cigar(cigar, cigarLength - 1 - c);
            char letter = ops_char[o.op];
            uint32_t length = o.len;

            switch (letter) {
                case '=':
                    aaIds += length;
                    // FALLTHROUGH
                case 'X':
                    // FALLTHROUGH
                case 'M':
                    queryPos += length;
                    targetPos += length;
                    backtrace.append(length, 'M');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'M');
                    break;
                case 'I':
                    queryPos += length;
                    backtrace.append(length, 'I');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'I');
                    break;
                case 'D':
                    targetPos += length;
                    backtrace.append(length, 'D');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'D');
                    break;
            }
            alnLength += length;
        }
    }

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.cigar = sAlnCigar;
    local_aln.cigarLen = cigarLength;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}


s_align
BlockAligner::bandedalign(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    s_align local_aln;
    // std::cout << "qIdx: " << qIdx << ", tIdx: " << tIdx << "\n";
    memset(sAlnCigar, 0, sizeof(uint32_t) * currentQuery->L);
    
    //Backward
    s_align local_aln_Backward;
    std::string backtrace_backward;
    if (qIdx > 0 && tIdx > 0) {
        local_aln_Backward = bandedalignBackward(currentTarget, qIdx, tIdx, backtrace_backward, xdrop);

    } else {
        local_aln_Backward.qStartPos1 = qIdx;
        local_aln_Backward.dbEndPos1 = tIdx;
        local_aln_Backward.score1 = 0;
        local_aln_Backward.identicalAACnt = 0;
    }

    //Forward
    std::string backtrace_forward;
    s_align local_aln_Forward = bandedalignForward(currentTarget, qIdx, tIdx, backtrace_forward, xdrop);
    if (qIdx > 0 && tIdx > 0) {
        // if backtrac_backward is empty-> print qIdx and tIdx
        if (backtrace_backward.empty()) {
            std::cout << "qIdx: " << qIdx << ", tIdx: " << tIdx << "\n";
        }
        if (backtrace_backward.back() != backtrace_forward.front()) {
            //remove one operation to avoid double counting the position at (qIdx, tIdx)
            std::cout << "backtrace_backward: " << backtrace_backward << "\n";
            std::cout << "backtrace_forward: " << backtrace_forward << "\n";
        }
    }
    //combine backtrace
    backtrace = backtrace_backward;;
    backtrace.append(backtrace_forward);

    local_aln.qStartPos1 = local_aln_Backward.qStartPos1;
    local_aln.dbStartPos1 = local_aln_Backward.dbStartPos1;
    local_aln.qEndPos1 = local_aln_Forward.qEndPos1;
    local_aln.dbEndPos1 = local_aln_Forward.dbEndPos1;
    local_aln.score1 = local_aln_Forward.score1 + local_aln_Backward.score1; // temporary
    local_aln.identicalAACnt = local_aln_Forward.identicalAACnt + local_aln_Backward.identicalAACnt -1 ; // temporary
    // sum of two cigars
    // uint32_t* concat = new uint32_t[size1 + size2];
    local_aln.cigar = local_aln_Forward.cigar; // temporary
    local_aln.cigarLen = local_aln_Backward.cigarLen + local_aln_Forward.cigarLen; // temporary

    float qcov = SmithWaterman::computeCov(local_aln.qStartPos1, local_aln.qEndPos1, currentQuery->L);
    float tcov = SmithWaterman::computeCov(local_aln.dbStartPos1, local_aln.dbEndPos1, currentTarget->L);
    double evalue = evaluer->computeEvalue(local_aln.score1, currentQuery->L);
    int bitScore = static_cast<int>(evaluer->computeBitScore(local_aln.score1) + 0.5);
    
    local_aln.qCov = qcov;
    local_aln.tCov = tcov;
    local_aln.evalue = evalue;
    local_aln.score1 = bitScore;

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    return local_aln;
}
