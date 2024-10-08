#include "nuclassembleUtil.h"
const double SMOOTHING_VALUE = 0.001;

#ifdef OPENMP
#include <omp.h>
#endif


void updateNuclAlignment(Matcher::result_t &tmpAlignment, DistanceCalculator::LocalAlignment &alignment,
                                const char *querySeq, size_t querySeqLen, const char *tSeq, size_t tSeqLen) {

    int qStartPos, qEndPos, dbStartPos, dbEndPos;
    int diag = alignment.diagonal;
    int dist = std::max(abs(diag), 0);

    if (diag >= 0) {
        qStartPos = alignment.startPos + dist;
        qEndPos = alignment.endPos + dist;
        dbStartPos = alignment.startPos;
        dbEndPos = alignment.endPos;
    } else {
        qStartPos = alignment.startPos;
        qEndPos = alignment.endPos;
        dbStartPos = alignment.startPos + dist;
        dbEndPos = alignment.endPos + dist;
    }

    int idCnt = 0;
    for(int i = qStartPos; i < qEndPos; i++){
        idCnt += (querySeq[i] == tSeq[dbStartPos+(i-qStartPos)]) ? 1 : 0;
    }
    float seqId =  static_cast<float>(idCnt) / (static_cast<float>(qEndPos) - static_cast<float>(qStartPos));

    tmpAlignment.seqId = seqId;
    tmpAlignment.qLen = querySeqLen;
    tmpAlignment.dbLen = tSeqLen;

    tmpAlignment.alnLength = alignment.diagonalLen;
    float scorePerCol = static_cast<float>(alignment.score ) / static_cast<float>(tmpAlignment.alnLength + 0.5);
    tmpAlignment.score = static_cast<int>(scorePerCol*100);

    tmpAlignment.qStartPos = qStartPos;
    tmpAlignment.qEndPos = qEndPos;
    tmpAlignment.dbStartPos = dbStartPos;
    tmpAlignment.dbEndPos = dbEndPos;

}

void getSeqErrorProf(diNucleotideProb & seqErrMatch, long double err)
{
    for (int orig = 0; orig < 4; orig++)
    {
        for (int obs = 0; obs < 4; obs++)
        {
            seqErrMatch.p[orig][obs] = (orig == obs) ? 1 - err : err/3;
        }
    }
    // for (int orig = 0; orig < 4; orig++)
    // {
    //     for (int obs = 0; obs < 4; obs++)
    //     {
    //         seqErrMis.p[orig][obs] = (( orig == 1 && obs == 3 ) || ( orig == 2 && obs == 0 )) ? 1 - err : err/3;
    //     }
    // }
}

char* getNuclRevFragment(char* fragment, size_t fragLen, NucleotideMatrix *nuclMatrix)
{
    char *fragmentRev = new char[fragLen];
    for (int pos = fragLen - 1; pos > -1; pos--) {
        int res = nuclMatrix->aa2num[static_cast<int>(fragment[pos])];
        char revRes = nuclMatrix->num2aa[nuclMatrix->reverseResidue(res)];
        fragmentRev[(fragLen - 1) - pos] = (revRes == 'X')? 'N' : revRes;
    }
    return fragmentRev;
}

float getRYSeqId(Matcher::result_t & res, char* querySeq,  char* targetSeq, std::unordered_map<char, int> & mapACGT)
{
    unsigned int distance = 0;
    unsigned int length = res.alnLength;

    // Check each position and count the differences 
    for (unsigned int t = 0; t < length; t++) {
        if ( mapACGT[querySeq[res.qStartPos + t]] != mapACGT[targetSeq[res.dbStartPos + t]] ) 
        {
            distance++;
        }
    }
    float rySeqid = static_cast<float>(res.alnLength - distance) / static_cast<float>(res.alnLength);
    return rySeqid;
}


void calcLikelihood(scorePerRes & scoredRes, char* querySeq, const char* targetSeq, std::vector<diNucleotideProb> & subDeamDiNuc, unsigned int maxAln, float randAlnPenal, diNucleotideProb & seqErrMatch, float excessPenal)
{
    // unsigned int countMatch = 0;
    // unsigned int countMismatch = 0;

    // Extract the query and target sequences
    std::string queryOverlap;
    for ( int i = scoredRes.r.qStartPos; i <= scoredRes.r.qEndPos; i++ )
    {
        queryOverlap += querySeq[i];
    }

    std::string targetOverlap;
    for ( int i = scoredRes.r.dbStartPos; i <= scoredRes.r.dbEndPos; i++ )
    {
        targetOverlap += targetSeq[i];
    }

    // Compare the sequences and print them if they're not identical
    /* if (queryOverlap != targetOverlap)
    {
        std::cerr << "isReverse: " << scoredRes.r.isRevToAlignment << std::endl;
        std::cerr << "Query Overlap: " << queryOverlap << std::endl;
        std::cerr << "targetSeq Overlap: " << targetOverlap << "\n" << std::endl;
    } */

    long double likMod = 0.0; // log(0)

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};


    std::vector<diNucleotideProb> subdeam_lookup(scoredRes.r.dbLen);
    // Create lookup vector
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[i];
    }
    for (size_t i = 5; i < scoredRes.r.dbLen - 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[5];
    }
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[scoredRes.r.dbLen - 5 + i] = subDeamDiNuc[ 6 + i];
    }

    //You either have a match or a mismatch
    for ( unsigned int t = 0; t < scoredRes.r.alnLength; t++ )
    {   
        double lik = 0;

        diNucleotideProb tProbs;
        tProbs = subdeam_lookup[scoredRes.r.dbStartPos + t];

        int qBase = nucleotideMap[queryOverlap[t]];
        int tBase = nucleotideMap[targetOverlap[t]];


        for (int target = 0; target<4; target++){
            double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qBase][target]);

            // seq error in obs. tBase:
            long double tBaseErr = 0;
            tBaseErr = seqErrMatch.p[target][tBase];

            //lik_one += std::exp(std::log(baseFreqs[query])+std::log(qBaseErr)+std::log(tBaseErr)+std::log(match_lik));
            lik += (tBaseErr * match_lik);
        }


    likMod += log(lik);

    } // else statment ending of the S or rate bases case senario


    // likelihood w/o penalty
    // double likAlnLen = likMod;

    // Penalizing short alignments
    unsigned int excess = maxAln - scoredRes.r.alnLength;
    likMod += ( excess * log(excessPenal) );

    if ( scoredRes.r.alnLength > maxAln )
    {
        std::cerr << "scoredRes.r.alnLength " <<  scoredRes.r.alnLength << std::endl;
        exit(1);
    }

    // Debugging
    if ( likMod > 0)
    {
        std::cerr << "Likelihood greater 0" << std::endl;
        exit(1);
    }

    // Calculate ratio (likMod / (likMod + random) )
    // In log space
    double randAln = maxAln * log(randAlnPenal);
    double ratioLog = 1.0/(1.0+exp(randAln-likMod));

    // Assign values to struct
    scoredRes.sLenNorm = likMod;
    scoredRes.sRatio = ratioLog;

} // END loop through each base


void calcLikelihoodConsensus(scorePerRes & scoredRes, std::string & consensus, const unsigned int querySeqLen, const char* targetSeq, std::vector<diNucleotideProb> & subDeamDiNuc, unsigned int & maxAln, float randAlnPenal, diNucleotideProb & seqErrMatch, float excessPenal)
{

    std::string targetOverlap;
    for ( unsigned int i = 0; i < scoredRes.r.dbLen; i++ )
    {
        targetOverlap += targetSeq[i];
    }

    long double likMod = 0.0; // log(0)

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};


    std::vector<diNucleotideProb> subdeam_lookup(scoredRes.r.dbLen);
    // Create lookup vector
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[i];
    }
    for (size_t i = 5; i < scoredRes.r.dbLen - 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[5];
    }
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[scoredRes.r.dbLen - 5 + i] = subDeamDiNuc[ 6 + i];
    }

    unsigned int alnCount = 0;

    unsigned int dbStartPos = scoredRes.r.dbStartPos;
    unsigned int dbEndPos = scoredRes.r.dbEndPos;
    unsigned int qStartPos = scoredRes.r.qStartPos;
    unsigned int qEndPos = scoredRes.r.qEndPos;

    const bool rightStart = dbStartPos == 0 && qEndPos == (querySeqLen - 1);
    const bool leftStart = qStartPos == 0 && dbEndPos == (scoredRes.r.dbLen - 1);

    if (leftStart){
        // padding
        unsigned int offset = scoredRes.r.dbLen - scoredRes.r.alnLength;
        unsigned int consensusStart = querySeqLen - offset;

        std::string leftPad(consensusStart, 'N');
        targetOverlap = leftPad + targetOverlap;

        unsigned int tIdx = 0;

        for (unsigned int i = 0; i < targetOverlap.size(); i++) {

            if ( !(targetOverlap[i] == 'N') ){
                tIdx++;
            }
            if (!(consensus[i] == 'N' || targetOverlap[i] == 'N')) {   
                alnCount++;
                double lik = 0;

                diNucleotideProb tProbs;
                tProbs = subdeam_lookup[tIdx - 1];

                int qBase = nucleotideMap[consensus[i]];
                int tBase = nucleotideMap[targetOverlap[i]];

                for (int target = 0; target<4; target++){
                    double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qBase][target]);

                    // seq error in obs. tBase:
                    long double tBaseErr = 0;
                    tBaseErr = seqErrMatch.p[target][tBase];

                    //lik_one += std::exp(std::log(baseFreqs[query])+std::log(qBaseErr)+std::log(tBaseErr)+std::log(match_lik));
                    lik += (tBaseErr * match_lik);
                }
                
                likMod += log(lik);
            }
        }

    } else if (rightStart){
        // padding
        unsigned int offset = scoredRes.r.dbLen - scoredRes.r.alnLength;
        unsigned int consensusStart = querySeqLen - offset;

        std::string rightPad(consensusStart, 'N');
        targetOverlap = targetOverlap + rightPad;

        unsigned int tIdx = 0;

        // Iterate from the end to the beginning
        for (unsigned int i = 0; i < targetOverlap.size(); i++) {
            
            unsigned int consensIdx = consensus.size() - targetOverlap.size() + i;
            if ( !(targetOverlap[i] == 'N') ){
                tIdx++;
            }
            if (!(consensus[consensIdx] == 'N' || targetOverlap[i] == 'N')) {
                alnCount++;
                double lik = 0;

                diNucleotideProb tProbs;
                tProbs = subdeam_lookup[tIdx - 1];

                int qBase = nucleotideMap[consensus[consensIdx]];
                int tBase = nucleotideMap[targetOverlap[i]];

                for (int target = 0; target<4; target++){
                    double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qBase][target]);

                    // seq error in obs. tBase:
                    long double tBaseErr = 0;
                    tBaseErr = seqErrMatch.p[target][tBase];

                    //lik_one += std::exp(std::log(baseFreqs[query])+std::log(qBaseErr)+std::log(tBaseErr)+std::log(match_lik));
                    lik += (tBaseErr * match_lik);

                }
                likMod += log(lik);
            }
        }
    }

    // likelihood w/o penalty
    // double likAlnLen = likMod;

    // Penalizing short alignments
    unsigned int excess = maxAln - alnCount;
    likMod += ( excess * log(excessPenal) );

    // Calculate ratio (likMod / (likMod + random) )
    // In log space
    double randAln = maxAln * log(randAlnPenal);
    double ratioLog = 1.0/(1.0+exp(randAln-likMod));

    // #pragma omp critical
    // {
    //     std::cerr << "c\t" << consensus << std::endl;
    //     std::cerr << "t\t" << targetOverlap << std::endl;
    //     std::cerr << "\n";
    //     std::cerr << "Target was extended?          " << scoredRes.r.wasExtended << std::endl;
    //     std::cerr << "                              log_space   exp_space" << std::endl;
    //     std::cerr << "target dbKey and len :        " << scoredRes.r.dbKey << " " << scoredRes.r.dbLen << std::endl;
    //     std::cerr << "query length                  " << scoredRes.r.qLen << std::endl;
    //     std::cerr << "likelihood no penalty :       " << likAlnLen << "   " << exp(likAlnLen) << std::endl;
    //     std::cerr << "likelihood with penalty:      " << likMod << "   " << exp(likMod) << std::endl;
    //     std::cerr << "random alignment:             " << randAln << "    " << exp(randAln) << std::endl;
    //     std::cerr << "                              ratio space" << std::endl;
    //     std::cerr << "ratio lik/(lik+random):       " << ratioLog  << std::endl;
    //     std::cerr << "alignment length:             " << scoredRes.r.alnLength << std::endl;
    //     std::cerr << "excess to longest:            " << excess << std::endl;
    //     //std::cerr << "ratio penal lik/(lik+random)  " << ratioLog2  << "  " << exp(ratioLog2) << std::endl;
    //     if ( alnCount > maxAln )
    //     {
    //         std::cerr << "alnCount " <<  alnCount << std::endl;
    //         exit(1);
    //     }

    //     // Debugging
    //     if ( likMod > 0)
    //     {
    //         std::cerr << "Likelihood greater 0" << std::endl;
    //         exit(1);
    //     }
    //     std::cerr << "\n";
    // }

    // Assign values to struct
    scoredRes.sLenNorm = likMod;
    scoredRes.sRatio = ratioLog;

}


void updateSeqIdConsensusReads(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, std::string & consensus, const unsigned int querySeqLen, unsigned int thread_idx, NucleotideMatrix *subMat, unsigned int & maxAlnLeft, unsigned int & maxAlnRight){
    // update sequence identity

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};

    std::unordered_map<char, int> ryMap = {
    {'A', 0},
    {'C', 1},
    {'G', 0},
    {'T', 1}};

    for (unsigned int idx = 0; idx < alignments.size(); idx++){
        // now retrieve the other candidate extensions and get their sequences

        unsigned int aln2updateId = sequenceDbr->getId(alignments[idx].dbKey);
        unsigned int aln2updateLen = sequenceDbr->getSeqLen(aln2updateId);

        char *aln2updateSequ = sequenceDbr->getData(aln2updateId, thread_idx);
        std::string aln2updateSeq;

        if (alignments[idx].isRevToAlignment) {
            // Convert the reversed fragment to std::string
            char *rightExtCandiSeqTmp = getNuclRevFragment(aln2updateSequ, aln2updateLen, (NucleotideMatrix *)subMat);
            aln2updateSeq = std::string(rightExtCandiSeqTmp, aln2updateLen);
            delete[] rightExtCandiSeqTmp;
        } else {
            // Directly convert the raw sequence to std::string
            aln2updateSeq = std::string(aln2updateSequ, aln2updateLen);
        }

        unsigned int dbStartPos = alignments[idx].dbStartPos;
        unsigned int dbEndPos = alignments[idx].dbEndPos;
        unsigned int qStartPos = alignments[idx].qStartPos;
        unsigned int qEndPos = alignments[idx].qEndPos;

        const bool rightStart = dbStartPos == 0 && qEndPos == (querySeqLen - 1);
        const bool leftStart = qStartPos == 0 && dbEndPos == (alignments[idx].dbLen - 1);

        int idCnt = 0;
        int idRyCnt = 0;
        unsigned int totalCnt = 0;
        if (leftStart){
            // padding
            unsigned int offset = alignments[idx].dbLen - alignments[idx].alnLength;
            unsigned int consensusStart = querySeqLen - offset;

            std::string leftPad(consensusStart, 'N');
            aln2updateSeq = leftPad + aln2updateSeq;

            for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
                if (!(consensus[i] == 'N' || aln2updateSeq[i] == 'N')) {
                    idCnt += (consensus[i] == aln2updateSeq[i]) ? 1 : 0;
                    idRyCnt += (ryMap[consensus[i]] == ryMap[aln2updateSeq[i]]) ? 1 : 0;
                    totalCnt += 1;
                }
            }
        } else if (rightStart){
            // padding
            unsigned int offset = alignments[idx].dbLen - alignments[idx].alnLength;
            unsigned int consensusStart = querySeqLen - offset;

            std::string rightPad(consensusStart, 'N');
            aln2updateSeq = aln2updateSeq + rightPad;

            for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
                unsigned int consensIdx = consensus.size() - aln2updateSeq.size() + i;
                if (!(consensus[consensIdx] == 'N' || aln2updateSeq[i] == 'N')) {
                    idCnt += (consensus[consensIdx] == aln2updateSeq[i]) ? 1 : 0;
                    idRyCnt += (ryMap[consensus[consensIdx]] == ryMap[aln2updateSeq[i]]) ? 1 : 0;
                    totalCnt += 1;
                }
            }
        }

        float seqId = alignments[idx].seqId;
        float rySeqId = alignments[idx].rySeqId;

        if (totalCnt != 0) {
            seqId = static_cast<float>(idCnt) / totalCnt;
            rySeqId = static_cast<float>(idRyCnt) / totalCnt;
        }

        if (leftStart && totalCnt > maxAlnLeft){
            maxAlnLeft = totalCnt;
        }
        else if (rightStart && totalCnt > maxAlnRight){
            maxAlnRight = totalCnt;
        }

        // #pragma omp critical
        // {
        //         std::cerr << "------------------------" << std::endl;
        //         std::cerr << "totalCnt " << totalCnt << std::endl;
        //         std::cerr << "idRyCnt " << idRyCnt << std::endl;
        //         std::cerr << "rySeqId " << rySeqId << std::endl;
        //         std::cerr << "idCnt " << idCnt << std::endl;
        //         std::cerr << "seqId " << seqId << std::endl;
        //         std::cerr << "query " << querySeq;
        //         std::cerr << "consensus " << consensus << std::endl;
        //         std::cerr << "aln2updateSeq " << aln2updateSeq << std::endl;
        //         std::cerr << "alignments[idx].seqId BEFORE " <<  alignments[idx].seqId << std::endl;
        //         std::cerr << "alignments[idx].rySeqId BEFORE " <<  alignments[idx].rySeqId << std::endl;
        //         std::cerr << "aln2updateId " <<  aln2updateId << std::endl;
        //         std::cerr << "queryKey " <<  queryKey << std::endl;
        //         std::cerr << "alignments[idx].dbKey " <<  alignments[idx].dbKey << std::endl;
        //         std::cerr << "alignments[idx].dbStart " <<  alignments[idx].dbStartPos << std::endl;
        //         std::cerr << "alignments[idx].dbEnd " <<  alignments[idx].dbEndPos << std::endl;
        //         std::cerr << "alignments[idx].qStart " <<  alignments[idx].qStartPos << std::endl;
        //         std::cerr << "alignments[idx].qEnd " <<  alignments[idx].qEndPos << std::endl;
        //         std::cerr << "alignments[idx].alnLength " <<  alignments[idx].alnLength << std::endl;
        //         std::cerr << "alignments[idx].dbLen " <<  alignments[idx].dbLen << std::endl;
        //         std::cerr << "queryseqlen " <<  querySeqLen << std::endl;
        //         std::cerr << "rightStart " <<  rightStart << std::endl;
        //         std::cerr << "leftStart " <<  leftStart << std::endl;
        // }

        alignments[idx].seqId = seqId;
        alignments[idx].rySeqId = rySeqId;
    }
}

// DEBUG TOOLS TO DELETE
void printDiNucleotideProb(const diNucleotideProb& matrix) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cerr << matrix.p[i][j] << "\t"; // Use tab for better spacing
        }
        std::cerr << std::endl;
    }
}
bool areDifferent(const diNucleotideProb& a, const diNucleotideProb& b) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (a.p[i][j] != b.p[i][j]) {
                return true;
            }
        }
    }
    return false;
}
void compareAndPrintIfDifferent(const diNucleotideProb& vec1, const diNucleotideProb& vec2) {

        if (areDifferent(vec1, vec2)) {
            std::cerr << "Matrix 1:" << std::endl;
            printDiNucleotideProb(vec1);
            std::cerr << "Matrix 2:" << std::endl;
            printDiNucleotideProb(vec2);
            std::exit(EXIT_FAILURE); // Exit the program
        }
}



// Function to calculate the consensus sequence
void calculateConsensus(std::string & consensus, const std::vector<std::vector<unsigned int>>& queryCov, unsigned int minCoverage) {
    // Initialize the consensus sequence with 'N's for the full length

    // Iterate over the range and build the consensus sequence
    for (size_t i = 0; i < queryCov.size(); ++i) {
        unsigned int total_coverage = queryCov[i][0] + queryCov[i][1] + queryCov[i][2] + queryCov[i][3];
        
        // Only consider positions with at least the minimum coverage
        if (total_coverage >= minCoverage) {
            unsigned int max_count = 0;
            char nucleotide = 'N'; // Default for positions without a clear consensus
            int count_of_max = 0; // To track how many times the max count is found

            // Determine the nucleotide with the highest count using switch-case
            for (int j = 0; j < 4; ++j) {
                if (queryCov[i][j] > max_count) {
                    max_count = queryCov[i][j];
                    nucleotide = "ACGT"[j];
                    count_of_max = 1; // Reset count of max occurrences
                } else if (queryCov[i][j] == max_count && max_count > 0) {
                    count_of_max++;
                }
            }

            // If the max count appears more than once, set nucleotide to 'N'
            if (count_of_max > 1) {
                nucleotide = 'N';
            }

            consensus[i] = nucleotide;
        }
    }
}


void consensusCaller(std::string & consensus, std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat){

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};

    std::unordered_map<char, int> ryMap = {
    {'A', 0},
    {'C', 1},
    {'G', 0},
    {'T', 1}};

    // std::cerr << "par.ancientUnsafe " << par.ancientUnsafe << std::endl;
    // exit(1);
    if ( par.ancientUnsafe == false ){
        for (unsigned int qPos = 0; qPos < querySeqLen; qPos++){
            unsigned int vecPos = querySeqLen + qPos;
            consensus[vecPos] = querySeq[qPos];  
        }
        return;
    }

    // 0. Initialize empty coverage vector
    // 1. Iterate through all alignments in the query
    // 2. Create the coverage vector for both ends
    // 3. Then find rule on how to decided between cases

    // Iterating through all alignments per query

    // Technically the query + the longest extension for both sides cannot excedd 3*querLen
    std::vector<std::vector<unsigned int>> queryCov(3*querySeqLen, std::vector<unsigned int> (4, 0));


    for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

        Matcher::result_t res = alignments[alnIdx];
        size_t dbKey = res.dbKey;
        unsigned int resId = sequenceDbr->getId(res.dbKey);

        // const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0 && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);


        // distinguish between right and left here

        unsigned int dbStartPos = res.dbStartPos;
        unsigned int dbEndPos = res.dbEndPos;
        unsigned int qStartPos = res.qStartPos;
        unsigned int qEndPos = res.qEndPos;

        if (resId == UINT_MAX) {
            Debug(Debug::ERROR) << "Could not find resId  " << res.dbKey
                                << " in database " << sequenceDbr->getDataFileName() << "\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int targetSeqLen = sequenceDbr->getSeqLen(resId);

        if ((rightStart || leftStart) && isNotIdentity){

            char *resSeq = sequenceDbr->getData(resId, thread_idx);
            unsigned int resSeqLen = sequenceDbr->getSeqLen(resId);
            std::string sequence;

            if ( res.isRevToAlignment ) {
                char *resSeqTmp = getNuclRevFragment(resSeq, res.dbLen, (NucleotideMatrix *) subMat);
                sequence = std::string(resSeqTmp, resSeqLen);
                delete[] resSeqTmp;
            } else {
                sequence = std::string(resSeq, resSeqLen);
            }

            if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                // right extension
                for (unsigned int pos = 0; pos < res.dbLen; pos ++){
                    unsigned int vecPos = querySeqLen + res.qStartPos + pos;
                    queryCov[vecPos][nucleotideMap[sequence[pos]]] += 1;
                }
                
            }
            else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                // left extension
                for (unsigned int pos = 0; pos < res.dbLen; pos ++){
                    unsigned int vecPos = querySeqLen - (res.dbLen - res.alnLength) + pos;
                    queryCov[vecPos][nucleotideMap[sequence[pos]]] += 1;
                }
            }
        }
    }

    calculateConsensus(consensus, queryCov, par.minCovSafe);

    // replace the query part of the consensus with the query since it is corrected and we trust it more than the consensus
    for (unsigned int qPos = 0; qPos < querySeqLen; qPos++){
        unsigned int vecPos = querySeqLen + qPos;
        consensus[vecPos] = querySeq[qPos];  
    }
    
    for (unsigned int qPos = 0; qPos < querySeqLen; qPos++){
        unsigned int vecPos = querySeqLen + qPos;
        queryCov[vecPos][nucleotideMap[querySeq[qPos]]] += 1;  
    }

    // #pragma omp critical
    // {   
    //     std::cerr << ">Consensus\n" << consensus << std::endl;
    //     std::cerr << ">Original\n" << querySeq;
    //     std::cerr << "QueryCov:\n";
    //     for (unsigned int i = 0; i < 3*querySeqLen; i++)
    //     {
    //         std::cerr << queryCov[i][0] << " ";
    //     }
    //     std::cerr << std::endl;
    //     for (unsigned int i = 0; i < 3*querySeqLen; i++)
    //     {
    //         std::cerr << queryCov[i][1] << " ";
    //     }
    //     std::cerr << std::endl;
    //     for (unsigned int i = 0; i < 3*querySeqLen; i++)
    //     {
    //         std::cerr << queryCov[i][2] << " ";
    //     }
    //     std::cerr << std::endl;
    //     for (unsigned int i = 0; i < 3*querySeqLen; i++)
    //     {
    //         std::cerr << queryCov[i][3] << " ";
    //     }
    //     std::cerr << "\n" << std::endl;
    // }
}

void updateSeqIdConsensus(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, std::string & consensus, const unsigned int querySeqLen, unsigned int thread_idx, NucleotideMatrix *subMat){
    // update sequence identity

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};

    std::unordered_map<char, int> ryMap = {
    {'A', 0},
    {'C', 1},
    {'G', 0},
    {'T', 1}};

    for (unsigned int idx = 0; idx < alignments.size(); idx++){
        // now retrieve the other candidate extensions and get their sequences

        unsigned int aln2updateId = sequenceDbr->getId(alignments[idx].dbKey);
        unsigned int aln2updateLen = sequenceDbr->getSeqLen(aln2updateId);

        char *aln2updateSequ = sequenceDbr->getData(aln2updateId, thread_idx);
        std::string aln2updateSeq;

        if (alignments[idx].isRevToAlignment) {
            // Convert the reversed fragment to std::string
            char *rightExtCandiSeqTmp = getNuclRevFragment(aln2updateSequ, aln2updateLen, (NucleotideMatrix *)subMat);
            aln2updateSeq = std::string(rightExtCandiSeqTmp, aln2updateLen);
            delete[] rightExtCandiSeqTmp;
        } else {
            // Directly convert the raw sequence to std::string
            aln2updateSeq = std::string(aln2updateSequ, aln2updateLen);
        }

        unsigned int dbStartPos = alignments[idx].dbStartPos;
        unsigned int dbEndPos = alignments[idx].dbEndPos;
        unsigned int qStartPos = alignments[idx].qStartPos;
        unsigned int qEndPos = alignments[idx].qEndPos;

        const bool rightStart = dbStartPos == 0 && qEndPos == (querySeqLen - 1);
        const bool leftStart = qStartPos == 0 && dbEndPos == (alignments[idx].dbLen - 1);

        int idCnt = 0;
        int idRyCnt = 0;
        int totalCnt = 0;
        if (leftStart){
            // padding
            unsigned int offset = alignments[idx].dbLen - alignments[idx].alnLength;
            unsigned int consensusStart = querySeqLen - offset;

            std::string leftPad(consensusStart, 'N');
            aln2updateSeq = leftPad + aln2updateSeq;

            for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
                if (!(consensus[i] == 'N' || aln2updateSeq[i] == 'N')) {
                    idCnt += (consensus[i] == aln2updateSeq[i]) ? 1 : 0;
                    idRyCnt += (ryMap[consensus[i]] == ryMap[aln2updateSeq[i]]) ? 1 : 0;
                    totalCnt += 1;
                }
            }
        } else if (rightStart){
            // padding
            unsigned int offset = alignments[idx].dbLen - alignments[idx].alnLength;
            unsigned int consensusStart = querySeqLen - offset;

            std::string rightPad(consensusStart, 'N');
            aln2updateSeq = aln2updateSeq + rightPad;

            // Iterate from the end to the beginning
            for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
                unsigned int consensIdx = consensus.size() - aln2updateSeq.size() + i;
                if (!(consensus[consensIdx] == 'N' || aln2updateSeq[i] == 'N')) {
                    idCnt += (consensus[consensIdx] == aln2updateSeq[i]) ? 1 : 0;
                    idRyCnt += (ryMap[consensus[consensIdx]] == ryMap[aln2updateSeq[i]]) ? 1 : 0;
                    totalCnt += 1;
                }
            }
        }

        float seqId = alignments[idx].seqId;
        float rySeqId = alignments[idx].rySeqId;

        if (totalCnt != 0) {
            seqId = static_cast<float>(idCnt) / totalCnt;
            rySeqId = static_cast<float>(idRyCnt) / totalCnt;
        }

        alignments[idx].seqId = seqId;
        alignments[idx].alnLengthCons = totalCnt;
        alignments[idx].rySeqId = rySeqId;
    }
    // #pragma omp critical
    // {
    //         std::cerr << "------------------------" << std::endl;
    //         // std::cerr << "totalCnt " << totalCnt << std::endl;
    //         // std::cerr << "idRyCnt " << idRyCnt << std::endl;
    //         // std::cerr << "rySeqId " << rySeqId << std::endl;
    //         // std::cerr << "idCnt " << idCnt << std::endl;
    //         // std::cerr << "seqId " << seqId << std::endl;
    //         std::cerr << "is unsafe " << par.ancientUnsafe << std::endl;
    //         std::cerr << "query " << querySeq;
    //         std::cerr << "consensus " << consensus << std::endl;
    //         // std::cerr << "aln2updateSeq " << aln2updateSeq << std::endl;
    //         // std::cerr << "alignments[idx].seqId BEFORE " <<  alignments[idx].seqId << std::endl;
    //         // std::cerr << "alignments[idx].rySeqId BEFORE " <<  alignments[idx].rySeqId << std::endl;
    //         // std::cerr << "aln2updateId " <<  aln2updateId << std::endl;
    //         std::cerr << "queryKey " <<  queryKey << std::endl;
    //         // std::cerr << "alignments[idx].dbKey " <<  alignments[idx].dbKey << std::endl;
    //         // std::cerr << "alignments[idx].dbStart " <<  alignments[idx].dbStartPos << std::endl;
    //         // std::cerr << "alignments[idx].dbEnd " <<  alignments[idx].dbEndPos << std::endl;
    //         // std::cerr << "alignments[idx].qStart " <<  alignments[idx].qStartPos << std::endl;
    //         // std::cerr << "alignments[idx].qEnd " <<  alignments[idx].qEndPos << std::endl;
    //         // std::cerr << "queryseqlen " <<  querySeqLen << std::endl;
    // }
}


void initDeamProbabilities(const std::string & deam5pfreqE,const std::string & deam3pfreqE, std::vector<substitutionRates> & sub5p, std::vector<substitutionRates> & sub3p, std::vector<diNucleotideProb> & allDeam, std::vector<diNucleotideProb> & revAllDeam){

    // std::vector<std::vector<diNucleotideProb>> & allDeam
    if ( deam3pfreqE == "3p.prof" && deam5pfreqE == "5p.prof"){
        for(unsigned int row=0;row<5;row++){
            substitutionRates noDamP;
            //fill.s[12]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            std::fill(std::begin(noDamP.s), std::end(noDamP.s), 0.0);
            sub5p.emplace_back(noDamP);
            sub3p.emplace_back(noDamP);
        }
    }
    else{
        readNucSubstitionRatesFreq(deam5pfreqE,sub5p);
        readNucSubstitionRatesFreq(deam3pfreqE,sub3p);
    }

    // printSubstitutionRatesVector(sub3p);
    // printSubstitutionRatesVector(sub5p);

    // Fill both ends with default deamination rates
    // Can later be adapted for flexible defualt values
    //double midSub5[16] = { 1, 0, 0, 0, 0, 0.95, 0, 0.05, 0, 0, 1, 0, 0, 0, 0, 1 };
    //double midSub3[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0.05, 0, 0.95, 0, 0, 0, 0, 1 };

    //float defaultBoth[16] = { 1, 0, 0, 0, 0, 0.99, 0, 0.01, 0.01, 0, 0.99, 0, 0, 0, 0, 1 };
    float defaultBoth[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };


    // Default matrix for middle part
    diNucleotideProb defaultDeam;
    // Iterate over all possible from-to combinations
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int idx = 4*i + j;
            defaultDeam.p[i][j] = defaultBoth[idx];
        }
    }

    // Create new "default-both-matrix" by taking over the last elements in the provided damage patterns
    if (!sub5p.empty()) {
        const auto& origStruct = sub5p.back(); // Directly access the last element
        // C->T: orig-struct.s[5]
        defaultDeam.p[1][3] = origStruct.s[5]; // C-T value
        defaultDeam.p[1][1] = 1 - origStruct.s[5]; // C-C value
    }

    if (!sub3p.empty()) {
        // Same here for G-A
        const auto& origStruct = sub3p.front(); // Directly access the last element
        // G->A: orig-struct.s[6]
        defaultDeam.p[2][0] = origStruct.s[6]; // G-A value
        defaultDeam.p[2][2] = 1 - origStruct.s[6]; // G-G value
    }

    std::unordered_map<int, double> defaultCT = {
        {5, defaultDeam.p[1][1]},
        {7, defaultDeam.p[1][3]}};
    std::unordered_map<int, double> defaultGA = {
        {8, defaultDeam.p[2][0]},
        {10, defaultDeam.p[2][2]}};

    // Initialize vector of diNucleotideProb structs
    std::vector<diNucleotideProb> subDeam;

    // Iterate over 5' original vector
    for (auto& origStruct : sub5p) {
        // Initialize new struct
        diNucleotideProb singleDeam;

        // Iterate over all possible from-to combinations
        int index5 = 0;
        for (int i = 0; i < 4; i++) {
            double sum = 0.0;
            for (int j = 0; j < 4; j++) {
                if (i == j)
                {
                    continue;
                }
                else {
                    singleDeam.p[i][j] = origStruct.s[index5];
                    sum += origStruct.s[index5];
                    index5++;
                }
            }
            singleDeam.p[i][i] = 1.0 - sum;
        }

        for (const auto& deam : defaultGA){
            int key = deam.first;
            int from = key/4;
            int to = key%4;
            singleDeam.p[from][to] = deam.second;
        }

        // Add new struct to vector
        subDeam.push_back(singleDeam);
    }


    // Iterate over 3' original vector
    for (auto& origStruct : sub3p) {
        // Initialize new struct
        diNucleotideProb singleDeam;

        // Iterate over all possible from-to combinations
        int index3 = 0;
        for (int i = 0; i < 4; i++) {
            double sum = 0.0;
            for (int j = 0; j < 4; j++) {
                if (i == j)
                {
                    continue;
                }
                else {
                    singleDeam.p[i][j] = origStruct.s[index3];
                    sum += origStruct.s[index3];
                    index3++;
                }
            }
            singleDeam.p[i][i] = 1.0 - sum;
        }

        for (const auto& deam : defaultCT){
            int key = deam.first;
            int from = key/4;
            int to = key%4;
            singleDeam.p[from][to] = deam.second;
        }

        // Add new struct to vector
        subDeam.push_back(singleDeam);
    }


    // Fill the first 5 positions with the first 5 elements of subDeam
    std::copy(subDeam.begin(), subDeam.begin() + 5, allDeam.begin());

    // Fill the last 5 positions with the last 5 elements of subDeam
    std::copy(subDeam.end() - 5, subDeam.end(), allDeam.end() - 5);

    // Fill the middle positions with defaultDeam
    std::fill(allDeam.begin() + 5, allDeam.end() - 5, defaultDeam);


    // create reverse that changes G-A values and C-T values:

    revAllDeam = allDeam;

    size_t deamLen = allDeam.size();

    for (size_t i = 0; i < deamLen; i++){
        diNucleotideProb end = allDeam[deamLen-1-i];

        revAllDeam[i].p[1][3] = end.p[2][0]; // revAllDeam.p[C][T] = end.p[G][A]
        revAllDeam[i].p[1][1] = end.p[2][2]; // revAllDeam.p[C][C] = end.p[G][G]

        revAllDeam[i].p[2][0] = end.p[1][3]; // revAllDeam.p[G][A] = begin.p[C][T]
        revAllDeam[i].p[2][2] = end.p[1][1]; // revAllDeam.p[G][G] = begin.p[C][C]

    }


    // std::cerr << "After initDeamProbabilities" << std::endl;

    // for ( unsigned int i = 0; i < deamLen; i++){
    //     std::cerr << "i: " << i << std::endl;
    //     std::cerr << "normal" << std::endl;
    //     for (int x = 0; x < 4; x++){
    //         for (int y = 0; y < 4; y++){
    //             std::cerr << allDeam[i].p[x][y] << " ";
    //         }
    //         std::cerr << std::endl;
    //     }
    //     std::cerr << std::endl;

    //     std::cerr << "rev: " << std::endl;
    //     for (int x = 0; x < 4; x++){
    //         for (int y = 0; y < 4; y++){
    //             std::cerr << revAllDeam[i].p[x][y] << " ";
    //         }
    //         std::cerr << std::endl;
    //     }
    //     std::cerr << std::endl;
    // } 

}//end initDeamProbabilities


// calculate number of matches/mismatches for adapted bayesian extension in ancientContigsResults.cpp
double deamMatches(Matcher::result_t res, unsigned int scoreAln, double matchLik){

    // Posterior:   p(match|C->T,d)
    // Likelihood:  p(C->T|match,d) = "from error profile"
    //              p(C->T|no match,d) = 5*(1-p(match,d))
    // Prior:       p(match) = 0.5 * ((score/ssingle_match)+0.9)/(L+1)) + 0.5 * alignment bias

    const double logAdjustmentConstant = std::log(1.4e-9);
    unsigned int maxLength = 1e5;

    // Helper lambda to compute log power adjustment
    auto logPower = [logAdjustmentConstant](unsigned int length) {
        return logAdjustmentConstant - 3.0 * std::log(length);
    };

    // Calculate log-adjustments at the boundaries
    double logMin = logPower(10); // Using max_length for min adjustment due to inversion
    double logMax = logPower(maxLength); // Using min_length for max adjustment due to inversion

    // Length calculation
    double logLength = logPower(std::min(res.alnLength, maxLength));
    double fractionLength = ( static_cast<double>(std::abs(logLength) - std::abs(logMax)) ) / static_cast<double>((std::abs(logMin) - std::abs(logMax)) );
    double priorAln = 1 - fractionLength;

    double pMatch = 0.5f * ((((static_cast<double>(scoreAln) + 3.0f * res.alnLength) / 5.0f) + 0.9f) / (res.alnLength + 1)) + 0.5f * priorAln;
    //double pMatch = ((((static_cast<float>(scoreAln) + 3.0f * res.alnLength) / 5.0f) + 0.9f) / (res.alnLength + 1)) * 0.9;

    double LikNoMatch = 1-pMatch;
    //double LikNoMatch = std::min(5*(1-pMatch), 0.25);
    double oddsRatio = LikNoMatch/matchLik;
    double odds = (1-pMatch)/pMatch;
    //double odds = (0.5)/pMatch;
    double posterior = 1 / (1 + oddsRatio * odds);
    // double posterior = pow((1 + (LikNoMatch/matchLik)*((1-pMatch)/pMatch)),-1);
 
    return posterior;
}


float ancientMatchCount(Matcher::result_t & res, std::string & consensus, const unsigned int querySeqLen, std::vector<diNucleotideProb> &subDeamDiNuc, DBReader<unsigned int>* sequenceDbr, unsigned int thread_idx, NucleotideMatrix *subMat, unsigned int & lefts, unsigned int & rights, unsigned int & others){

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};

    std::unordered_map<char, int> ryMap = {
    {'A', 0},
    {'C', 1},
    {'G', 0},
    {'T', 1}};


    // Final match calculation
    // m = 1/2*(L+(score/ssingle_match))+sum(p(C->T|match,d)-ssingle_mis/ssingle_match)+sum(p(G->A|match,d)-ssingle_mis/ssingle_match))
    // const bool hasOverlap = res.dbLen > res.alnLength; 

    float mCT = 0;
    float mGA = 0;

    unsigned int mmCons = (1 - res.seqId) * res.alnLengthCons + 0.5;
    unsigned int mCons = res.alnLengthCons - mmCons;
    unsigned int scoreAln = mCons * 2 + mmCons * (-3); // TODO: use match score from nucleotide.out


    unsigned int aln2updateId = sequenceDbr->getId(res.dbKey);
    unsigned int aln2updateLen = sequenceDbr->getSeqLen(aln2updateId);

    char *aln2updateSequ = sequenceDbr->getData(aln2updateId, thread_idx);
    std::string aln2updateSeq;

    if (res.isRevToAlignment) {
        // Convert the reversed fragment to std::string
        char *rightExtCandiSeqTmp = getNuclRevFragment(aln2updateSequ, aln2updateLen, (NucleotideMatrix *)subMat);
        aln2updateSeq = std::string(rightExtCandiSeqTmp, aln2updateLen);
        delete[] rightExtCandiSeqTmp;
    } else {
        // Directly convert the raw sequence to std::string
        aln2updateSeq = std::string(aln2updateSequ, aln2updateLen);
    }


    unsigned int dbStartPos = res.dbStartPos;
    unsigned int dbEndPos = res.dbEndPos;
    unsigned int qStartPos = res.qStartPos;
    unsigned int qEndPos = res.qEndPos;

    const bool rightStart = dbStartPos == 0 && qEndPos == (querySeqLen - 1);
    const bool leftStart = qStartPos == 0 && dbEndPos == (res.dbLen - 1);

    // count left and right start and others.

    // std::cerr << "Before left Start" << std::endl;
    if (leftStart){
        // padding

        lefts++;

        unsigned int offset = res.dbLen - res.alnLength;
        unsigned int consensusStart = querySeqLen - offset;

        std::string leftPad(consensusStart, 'N');
        aln2updateSeq = leftPad + aln2updateSeq;

        for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
            if (!(consensus[i] == 'N' || aln2updateSeq[i] == 'N')) {

                int qBase = nucleotideMap[consensus[i]];
                int tBase = nucleotideMap[aln2updateSeq[i]];
                
                // check for C-to-T and G-to-A:
                bool dimerCT = ((4*qBase + tBase) == 7 ); // = 4*1 + 3
                bool dimerGA = ((4*qBase + tBase) == 8 ); // = 4*2 + 0

                double matchLik = 0;
                matchLik = subDeamDiNuc[5].p[qBase][tBase];

                if ( dimerCT && matchLik > 0 ){
                    mCT += deamMatches(res, scoreAln, matchLik);
                }
                else if ( dimerGA && matchLik > 0 ){
                    //std::cerr << "matchLik GA\t" << matchLik << std::endl;
                    mGA += deamMatches(res, scoreAln, matchLik);
                }
            }
        }
    // std::cerr << "End left Start" << std::endl;
    // std::cerr << "Before right Start" << std::endl;
    } else if (rightStart){

        rights++;
        // padding
        unsigned int offset = res.dbLen - res.alnLength;
        unsigned int consensusStart = querySeqLen - offset;

        std::string rightPad(consensusStart, 'N');
        aln2updateSeq = aln2updateSeq + rightPad;

        // Iterate from the end to the beginning
        for (unsigned int i = 0; i < aln2updateSeq.size(); i++) {
            unsigned int consensIdx = consensus.size() - aln2updateSeq.size() + i;
            if (!(consensus[consensIdx] == 'N' || aln2updateSeq[i] == 'N')) {

                int qBase = nucleotideMap[consensus[consensIdx]];
                int tBase = nucleotideMap[aln2updateSeq[i]];
                
                // check for C-to-T and G-to-A:
                bool dimerCT = ((4*qBase + tBase) == 7 ); // = 4*1 + 3
                bool dimerGA = ((4*qBase + tBase) == 8 ); // = 4*2 + 0

                double matchLik = 0;
                matchLik = subDeamDiNuc[5].p[qBase][tBase];

                if ( dimerCT && matchLik > 0 ){
                    mCT += deamMatches(res, scoreAln, matchLik);
                }
                else if ( dimerGA && matchLik > 0 ){
                    //std::cerr << "matchLik GA\t" << matchLik << std::endl;
                    mGA += deamMatches(res, scoreAln, matchLik);
                }
            }
        }
    }
    else{
        others++;
    }
    // std::cerr << "End right Start" << std::endl;

    float matches = ((static_cast<float>(scoreAln) + 3.0f*res.alnLengthCons) / 5.0f) + mCT + mGA; // TODO: use match score from nucleotide.out
    return matches;
}