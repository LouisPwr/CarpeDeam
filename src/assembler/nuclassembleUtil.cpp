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


float ancientMatchCount(Matcher::result_t res, char* querySeq, char* targetSeq, std::vector<diNucleotideProb> &subDeamDiNuc, std::unordered_map<char, int> &nucleotideMap){

    // Final match calculation
    // m = 1/2*(L+(score/ssingle_match))+sum(p(C->T|match,d)-ssingle_mis/ssingle_match)+sum(p(G->A|match,d)-ssingle_mis/ssingle_match))
    // const bool hasOverlap = res.dbLen > res.alnLength; 

    float mCT = 0;
    float mGA = 0;

    unsigned int mmCons = (1 - res.seqId) * res.alnLength + 0.5;
    unsigned int mCons = res.alnLength - mmCons;
    unsigned int scoreAln = mCons * 2 + mmCons * (-3); // TODO: use match score from nucleotide.out

    for (unsigned int pos = 0; pos < res.alnLength; pos++){
        
        int qBase = nucleotideMap[querySeq[res.qStartPos + pos]];
        int tBase = nucleotideMap[targetSeq[res.dbStartPos + pos]];
        
        // check for C-to-T and G-to-A:
        bool dimerCT = ((4*qBase + tBase) == 7 ); // = 4*1 + 3
        bool dimerGA = ((4*qBase + tBase) == 8 ); // = 4*2 + 0

        double matchLik = 0;
        //const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        //const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);
        matchLik = subDeamDiNuc[5].p[qBase][tBase];

        if ( dimerCT && matchLik > 0 ){
            mCT += deamMatches(res, scoreAln, matchLik);
        }
        else if ( dimerGA && matchLik > 0 ){
            //std::cerr << "matchLik GA\t" << matchLik << std::endl;
            mGA += deamMatches(res, scoreAln, matchLik);
        }
    }
    float matches = ((static_cast<float>(scoreAln) + 3.0f*res.alnLength) / 5.0f) + mCT + mGA; // TODO: use match score from nucleotide.out
    return matches;
}


float ancientMatchCountContig(Matcher::result_t res, char* querySeq, char* targetSeq, std::vector<diNucleotideProb> &subDeamDiNuc, std::unordered_map<char, int> &nucleotideMap){

    // Final match calculation
    // m = 1/2*(L+(score/ssingle_match))+sum(p(C->T|match,d)-ssingle_mis/ssingle_match)+sum(p(G->A|match,d)-ssingle_mis/ssingle_match))
    // const bool hasOverlap = res.dbLen > res.alnLength; 

    float mCT = 0;
    float mGA = 0;
    unsigned int mmCons = (1 - res.seqId) * res.alnLength + 0.5;
    unsigned int mCons = res.alnLength - mmCons;
    unsigned int scoreAln = mCons * 2 + mmCons * (-3); // TTODO: use match score from nucleotide.out
    for (unsigned int pos = 0; pos < res.alnLength; pos++){
        
        int qBase = nucleotideMap[querySeq[res.qStartPos + pos]];
        int tBase = nucleotideMap[targetSeq[res.dbStartPos + pos]];

        // check for C-to-T and G-to-A:
        bool dimerCT = ((4*qBase + tBase) == 7 ); // = 4*1 + 3
        bool dimerGA = ((4*qBase + tBase) == 8 ); // = 4*2 + 0

        double matchLik = subDeamDiNuc[5].p[qBase][tBase];
        
        if ( dimerCT && matchLik > 0 ){
            //std::cerr << "matchLik CT\t" << matchLik << std::endl;
            mCT += deamMatches(res, scoreAln, matchLik);
            //std::cerr << "mCT " << mCT << std::endl;
        }
        else if ( dimerGA && matchLik > 0 ){
            //std::cerr << "matchLik GA\t" << matchLik << std::endl;
            mGA += deamMatches(res, scoreAln, matchLik);
        }
    }
    float matches = 0.5*(res.alnLength + static_cast<float>(scoreAln)/static_cast<float>(2) + mCT + mGA); // TODO: use match score from nucleotide.out
    return matches;
}


float getSubSeqId(Matcher::result_t &aln, const char *querySeq, const char *tSeq) {

    int qStart = aln.qStartPos + 5;
    int qEnd = aln.qEndPos - 5;
    int tStart = aln.dbStartPos + 5;
    int idCnt = 0;
    
    for(int i = qStart; i < qEnd; i++){
        idCnt += (querySeq[i] == tSeq[tStart+(i-qStart)]) ? 1 : 0;
    }
    float seqId =  static_cast<float>(idCnt) / (static_cast<float>(qEnd) - static_cast<float>(qStart));
    seqId = std::round(seqId * 1000.0f) / 1000.0f;
    return seqId;
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

    // Here we start calculating the likelihood for a >>> S I N G L E <<< alignment
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


    // for (unsigned int t = 0; t < scoredRes.r.alnLength; t++) {
    //     long double lik = 0;
        
    //     diNucleotideProb tProbs = subdeam_lookup[scoredRes.r.dbStartPos + t];
        
    //     int qBase = nucleotideMap[queryOverlap[t]];
    //     int tBase = nucleotideMap[targetOverlap[t]];
        
            
    //     for (int target = 0; target < 4; target++) {
    //         double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qBase][target]);
            
    //         // seq error in obs. tBase:
    //         long double tBaseErr = seqErrMatch.p[target][tBase];

    //         lik += (tBaseErr * match_lik);
    //     }
        
    //     likMod += std::log(lik);
    // }

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
    double likAlnLen = likMod;

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
    //double randAln = countMismatch * log(1.0/2.0) + countMatch * log(1.0/2.0);
    double randAln = maxAln * log(randAlnPenal);
    //scoredRes.r.alnLength * log(3.0/4.0);
    // Now in log space
    //double ratio = log(exp(likAlnLen)) - log( exp(likAlnLen) + exp(randAln) );
    double ratioLog = 1.0/(1.0+exp(randAln-likMod));
    //double ratioLog2 = 1.0/(1.0+exp(likMod-randAln));

    if ( false )
    //if ( randAln > likMod )
{
        std::cerr << "q" << queryOverlap << std::endl;
        std::cerr << "t" << targetOverlap << std::endl;
        std::cerr << "\n";
        std::cerr << "Target was extended?          " << scoredRes.r.wasExtended << std::endl;
        std::cerr << "                              log_space   exp_space" << std::endl;
        std::cerr << "target dbKey and len :        " << scoredRes.r.dbKey << " " << scoredRes.r.dbLen << std::endl;
        std::cerr << "query length                  " << scoredRes.r.qLen << std::endl;
        std::cerr << "likelihood no penalty :       " << likAlnLen << "   " << exp(likAlnLen) << std::endl;
        std::cerr << "likelihood with penalty:      " << likMod << "   " << exp(likMod) << std::endl;
        std::cerr << "random alignment:             " << randAln << "    " << exp(randAln) << std::endl;
        std::cerr << "                              ratio space" << std::endl;
        std::cerr << "ratio lik/(lik+random):       " << ratioLog  << std::endl;
        std::cerr << "alignment length:             " << scoredRes.r.alnLength << std::endl;
        std::cerr << "excess to longest:            " << excess << std::endl;
        //std::cerr << "ratio penal lik/(lik+random)  " << ratioLog2  << "  " << exp(ratioLog2) << std::endl;
        std::cerr << "\n";
}    

    // Assign values to struct
    scoredRes.sLenNorm = likMod;
    scoredRes.sRatio = ratioLog;

} // END loop through each base


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

// // Function to calculate the consensus sequence
// std::string calculateConsensus(const std::vector<std::vector<unsigned int>>& queryCov, unsigned int minCoverage = 3) {
//     // Initialize the consensus sequence with 'N's for the full length
//     std::string consensus(queryCov.size(), 'N');

//     // Iterate over the range and build the consensus sequence
//     for (size_t i = 0; i < queryCov.size(); ++i) {
//         unsigned int total_coverage = queryCov[i][0] + queryCov[i][1] + queryCov[i][2] + queryCov[i][3];
        
//         // Only consider positions with at least the minimum coverage
//         if (total_coverage >= minCoverage) {
//             unsigned int max_count = 0;
//             char nucleotide = 'N'; // Default for positions without a clear consensus
            
//             // Determine the nucleotide with the highest count using switch-case
//             for (int j = 0; j < 4; ++j) {
//                 if (queryCov[i][j] > max_count) {
//                     max_count = queryCov[i][j];
//                     switch (j) {
//                         case 0: nucleotide = 'A'; break;
//                         case 1: nucleotide = 'C'; break;
//                         case 2: nucleotide = 'G'; break;
//                         case 3: nucleotide = 'T'; break;
//                     }
//                 }
//             }
//             consensus[i] = nucleotide;
//         }
//     }
//     return consensus;
// }


// Function to calculate the consensus sequence
std::string calculateConsensus(const std::vector<std::vector<unsigned int>>& queryCov, unsigned int minCoverage = 2) {
    // Initialize the consensus sequence with 'N's for the full length
    std::string consensus(queryCov.size(), 'N');

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

    return consensus;
}


std::string consensusCaller(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat){

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

            // #pragma omp critical
            // {
            // std::cerr << "query " << querySeq << std::endl;
            // std::cerr << "sequence " << sequence << std::endl;
            // }

            if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                // right extension
                for (unsigned int pos = 0; pos < res.dbLen; pos ++){
                    unsigned int vecPos = querySeqLen - 1 + res.qStartPos + pos;
                    queryCov[vecPos][nucleotideMap[sequence[pos]]] += 1;
                }
                
            }
            else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                // left extension
                for (unsigned int pos = 0; pos < res.dbLen; pos ++){
                    unsigned int vecPos = querySeqLen - 1 - (res.dbLen - res.alnLength) + pos;
                    queryCov[vecPos][nucleotideMap[sequence[pos]]] += 1;
                }
            }
        }
    }

    std::string consensus = calculateConsensus(queryCov);

    // replace the query part of the consensus with the query since it is corrected and we trust it more than the consensus
    for (unsigned int qPos = 0; qPos < querySeqLen; qPos++){
        unsigned int vecPos = querySeqLen - 1 + qPos;
        consensus[vecPos] = querySeq[qPos];  
    }
    
    for (unsigned int qPos = 0; qPos < querySeqLen; qPos++){
        unsigned int vecPos = querySeqLen - 1 + qPos;
        queryCov[vecPos][nucleotideMap[querySeq[qPos]]] += 1;  
    }

    return consensus;

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

void updateSeqIdConsensus(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, std::string & consensus, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat){
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
            unsigned int consensusStart = querySeqLen - offset - 1;

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
            unsigned int consensusStart = querySeqLen - offset + 1;

            std::string rightPad(consensusStart, 'N');
            aln2updateSeq = aln2updateSeq + rightPad;

            // Iterate from the end to the beginning
            for (unsigned int i = aln2updateSeq.size() - 1, j = consensus.size() - 1; i == 0 && j == 0; --i, --j) {
                if (!(consensus[j] == 'N' || aln2updateSeq[i] == 'N')) {
                    idCnt += (consensus[j] == aln2updateSeq[i]) ? 1 : 0;
                    idRyCnt += (ryMap[consensus[j]] == ryMap[aln2updateSeq[i]]) ? 1 : 0;
                    totalCnt += 1;
                }
                // Ensure that we do not underflow the size_t type
                if (i == 0 || j == 0) break;
            }
        }

        float seqId = alignments[idx].seqId;
        float rySeqId = alignments[idx].rySeqId;

        if (totalCnt != 0) {
            seqId = static_cast<float>(idCnt) / totalCnt;
            rySeqId = static_cast<float>(idRyCnt) / totalCnt;
        }

        // #pragma omp critical
        // {
        //     if ( queryKey == 4631784 ){
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
        //         std::cerr << "queryseqlen " <<  querySeqLen << std::endl;
        //     }
        // }

        alignments[idx].seqId = seqId;
        alignments[idx].rySeqId = rySeqId;
    }
}




std::pair<bool, bool> mgeFinderContigs(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat){

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

// MGE marker
    bool mgeFoundRight = false;
    bool mgeFoundLeft = false;

    std::vector<Matcher::result_t> mgeCandiRight;
    std::vector<Matcher::result_t> mgeCandiLeft;

    float rymerThresh = par.correctionThresholdRySeqId;

    // define thresholds to include extension candidates for the mge-identifier
    // float alnSeqIdThr = par.seqIdThr;
    //float alnSeqIdThr = 0.999;
    //unsigned int extLen = 150;

    //set threshold when to declare an mge
    // float mgeSeqId = 0.9;
    // float mgeRySeqId = 0.99;

    //set number of bases in overlap to comapare after the alignment ends
    //unsigned int numBaseCompare = 50;

    for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

        Matcher::result_t res = alignments[alnIdx];
        size_t dbKey = res.dbKey;
        unsigned int resId = sequenceDbr->getId(res.dbKey);
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0 && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);


        // distinguish between right and left here

        unsigned int dbStartPos = res.dbStartPos;
        unsigned int dbEndPos = res.dbEndPos;
        unsigned int qStartPos = res.qStartPos;
        unsigned int qEndPos = res.qEndPos;
        unsigned int targetId = sequenceDbr->getId(res.dbKey);
        if (targetId == UINT_MAX) {
            Debug(Debug::ERROR) << "Could not find targetId  " << res.dbKey
                                << " in database " << sequenceDbr->getDataFileName() << "\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId);

        // TODO: Add min overlap len and max target seq len
        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){

            char *resSeq = sequenceDbr->getData(resId, thread_idx); 
            std::string resSeqStr;
            float ryId;

            if ( res.alnLength <= 100){
                rymerThresh = (static_cast<float>(res.alnLength) - 1) / static_cast<float>(res.alnLength);
                rymerThresh = std::floor(rymerThresh * 1000) / 1000;
            }

            if ( res.isRevToAlignment ) {
                char *resSeqTmp = getNuclRevFragment(resSeq, res.dbLen, (NucleotideMatrix *) subMat);
                ryId = getRYSeqId(res, querySeq, resSeqTmp, ryMap);
                res.rySeqId = ryId;
                delete[] resSeqTmp;
            }
            else {
                ryId = getRYSeqId(res, querySeq, resSeq, ryMap);
                res.rySeqId = ryId;
            }

            if ( res.seqId >= 0.99 ){
                // count only the left and right extensions which extensions should be compared
                //if (res.seqId > alnSeqIdThr && (res.dbLen - res.alnLength <= extLen)){
                //if (res.seqId >= alnSeqIdThr){
                if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                    // right extension
                    mgeCandiRight.push_back(res);
                }
                else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                    // left extension
                    mgeCandiLeft.push_back(res);
                }
                //}
            }
        }
    }


#ifdef DEBUGCORR

        int leadingDashes = 30;
        if ( queryKey == 4631784 ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string headers = "sequences\toverlap_similarity\toverlap_rySimilarity\tsimilarity_to_longest_extension\trySimilarity_to_longest_extension\texpected_similarity";

                std::string qeueryDebug;
                // Construct qeueryDebug from querySeq
                for (int i = 0; i < querySeqLen; i++) {
                    qeueryDebug += querySeq[i];
                }

                // Create a string of 100 "-" characters
                std::string dashes(leadingDashes, '*');

                // Add dashes before and after qeueryDebug
                qeueryDebug = dashes + qeueryDebug + dashes;

                if (outputFileOver.is_open()) {
                    outputFileOver << headers << "\n";
                    outputFileOver << qeueryDebug << "\n";
                }
            }
        }
#endif

    // Now we iterate through the left and right extensions
    // Get the longest extension and align all remaining candidates
    // right extensions
    if (mgeCandiRight.size() > 1){
        
        // Get extension candidate with the longest extension 
        Matcher::result_t rightLongestExt = mgeCandiRight.front();
        for (unsigned int rightCandis = 0; rightCandis < mgeCandiRight.size(); rightCandis++){
            rightLongestExt = (mgeCandiRight[rightCandis].dbLen - mgeCandiRight[rightCandis].alnLength >= rightLongestExt.dbLen - rightLongestExt.alnLength) ? mgeCandiRight[rightCandis] : rightLongestExt;
        }

        // now retrieve the sequence of the longest extension in the correct orientation 
        unsigned int rightLongestExtId = sequenceDbr->getId(rightLongestExt.dbKey);
        unsigned int rightLongestExtLen = sequenceDbr->getSeqLen(rightLongestExtId);

        char *rightLongestExtSequ = sequenceDbr->getData(rightLongestExtId, thread_idx);
        std::string rightLongestExtSeq;

        if ( rightLongestExt.isRevToAlignment ) {
            char *rightLongestExtSeqTmp = getNuclRevFragment(rightLongestExtSequ, rightLongestExtLen, (NucleotideMatrix *) subMat);
            rightLongestExtSeq = std::string(rightLongestExtSeqTmp, rightLongestExtLen);
            delete[] rightLongestExtSeqTmp;
        }
        else {
            rightLongestExtSeq = std::string(rightLongestExtSequ, rightLongestExtLen);
        }
        
#ifdef DEBUGCORR
            if ( queryKey == 4631784 ){
            // DEBUG CORRECTION
                #pragma omp critical
                {
                    std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '+');

                    for (int i = leadingDashes + (querySeqLen - rightLongestExt.alnLength), k = 0; k < rightLongestExt.dbLen; i++, k++) {
                            targetAligned[i] = rightLongestExtSeq[k];
                    }      

                    if (outputFileOver.is_open()) {
                        outputFileOver << targetAligned << "\n";
                    }
                }
            }
#endif

        // Iterate through all candidates and compare to the longest
        for (unsigned int rightRes = 0; rightRes < mgeCandiRight.size(); rightRes++){
            // now retrieve the other candidate extensions and get their sequences

            Matcher::result_t righty = mgeCandiRight[rightRes];

            unsigned int rightExtCandiId = sequenceDbr->getId(righty.dbKey);
            unsigned int rightExtCandiLen = sequenceDbr->getSeqLen(rightExtCandiId);

            if ( rightExtCandiId == rightLongestExtId ){
                continue;
            }

            char *rightExtCandiSequ = sequenceDbr->getData(rightExtCandiId, thread_idx);
            std::string rightExtCandiSeq;

            //bool isBadCandi = false;
            //std::vector<long double> expSim;
            if (righty.isRevToAlignment) {
                // Convert the reversed fragment to std::string
                char *rightExtCandiSeqTmp = getNuclRevFragment(rightExtCandiSequ, rightExtCandiLen, (NucleotideMatrix *)subMat);
                rightExtCandiSeq = std::string(rightExtCandiSeqTmp, rightExtCandiLen);
                delete[] rightExtCandiSeqTmp;
                // Here we calculate the epected similarity of the extension
                // expSim = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNucRev, seqErrMatch, true);
            } else {
                // Directly convert the raw sequence to std::string
                rightExtCandiSeq = std::string(rightExtCandiSequ, rightExtCandiLen);
                // Here we calculate the epected similarity of the extension
                // expSim = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNuc, seqErrMatch, true);
            }

            // calculate the sequence identity
            // unsigned int othersExtLen = ( righty.dbLen - righty.alnLength < numBaseCompare ) ? righty.dbLen - righty.alnLength : numBaseCompare;
            unsigned int othersExtLen = righty.dbLen - righty.alnLength;
            int idCnt = 0;
            int idRyCnt = 0;
            for (unsigned int rightPos = 0; rightPos < othersExtLen; ++rightPos) {
                idCnt += (rightLongestExtSeq[rightLongestExt.dbEndPos + rightPos] == rightExtCandiSeq[righty.dbEndPos + rightPos]) ? 1 : 0;
                idRyCnt += (ryMap[rightLongestExtSeq[rightLongestExt.dbEndPos + rightPos]] == ryMap[rightExtCandiSeq[righty.dbEndPos + rightPos]]) ? 1 : 0;
            }

            //float seqId = static_cast<float>(idCnt) / othersExtLen;
            float rySeqId = static_cast<float>(idRyCnt) / othersExtLen;

                // #pragma omp critical
                // {
                //     std::cerr << "MGEFINDER RIGHT" << std::endl;
                //     std::cerr << "othersExtLen " << othersExtLen << std::endl;
                //     std::cerr << "idRyCnt " << idRyCnt << std::endl;
                //     std::cerr << "rySeqId " << rySeqId << std::endl;
                //     std::cerr << "idCnt " << idCnt << std::endl;
                //     std::cerr << "seqId " << seqId << std::endl;
                //     std::cerr << "query " << querySeq << std::endl;
                //     std::cerr << "rightLongestExtSeq " << rightLongestExtSeq << std::endl;
                //     std::cerr << "rightExtCandiSeq " << rightExtCandiSeq << std::endl;
                // }

            // check for number of C-T and G-A mismatches:
            unsigned int ryMM = idRyCnt - idCnt;
            unsigned int ryMMThr = othersExtLen * 0.1; //allowing for 2% mismatch rate

            // if ( rySeqId < rymerThresh ){
            // //if ( rySeqId < rymerThresh && ( (othersExtLen <= 10 && ryMM <= 2) || ryMM <= ryMMThr )){
            // //if ( rySeqId < rymerThresh || ( (othersExtLen < 10 && ryMM > 2) || seqId < 0.95 )){
            // //if ( rySeqId < rymerThresh || ( (othersExtLen < 10 && ryMM > 2) || seqId < 0.95 )){
            //     mgeFoundRight = true;
            //     break;
            // } 

            if (othersExtLen <= 10 && ryMM > 1) {
                mgeFoundRight = true;
                break;
            } else if (rySeqId < rymerThresh) {
                mgeFoundRight = true;
                break;
            } else if ( ryMM > ryMMThr ){
                mgeFoundRight = true;
                break;  
            }

#ifdef DEBUGCORR
        if ( queryKey == 4631784 ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '-');
            
                for (int i = leadingDashes + (querySeqLen - righty.alnLength), k = 0; k < righty.dbLen; i++, k++) {
                    targetAligned[i] = rightExtCandiSeq[k];
                }

                // Append the loop output to outputLine
                targetAligned += "\t";
                targetAligned += std::to_string(righty.seqId);
                targetAligned += "\t";
                // targetAligned += std::to_string(insideSeqId);
                // targetAligned += "\t";
                targetAligned += std::to_string(righty.rySeqId);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[0]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[1]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[2]);
                targetAligned += "\n";

                if (outputFileOver.is_open()) {
                    outputFileOver << targetAligned;
                }
            }
        }
#endif
        }
    }

    // left extensions
    if (mgeCandiLeft.size() > 1){
        Matcher::result_t leftLongestExt = mgeCandiLeft.front();

        for (unsigned int leftCandis = 0; leftCandis < mgeCandiLeft.size(); leftCandis++){
            leftLongestExt = (mgeCandiLeft[leftCandis].dbLen - mgeCandiLeft[leftCandis].alnLength >= leftLongestExt.dbLen - leftLongestExt.alnLength) ? mgeCandiLeft[leftCandis] : leftLongestExt;
        }

        // now retrieve the sequence of the longest extension in the correct orientation 
        unsigned int leftLongestExtId = sequenceDbr->getId(leftLongestExt.dbKey);
        unsigned int leftLongestExtLen = sequenceDbr->getSeqLen(leftLongestExtId);

        char *leftLongestExtSequ = sequenceDbr->getData(leftLongestExtId, thread_idx);
        std::string leftLongestExtSeq;

        if ( leftLongestExt.isRevToAlignment ) {
            char *leftLongestExtSeqTmp = getNuclRevFragment(leftLongestExtSequ, leftLongestExtLen, (NucleotideMatrix *) subMat);
            leftLongestExtSeq = std::string(leftLongestExtSeqTmp, leftLongestExtLen);
            delete[] leftLongestExtSeqTmp;
        }
        else {
            leftLongestExtSeq = std::string(leftLongestExtSequ, leftLongestExtLen);
        }

#ifdef DEBUGCORR
            if ( queryKey == 4631784 ){
            // DEBUG CORRECTION
                #pragma omp critical
                {
                    std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '+');

                    for (int i = leadingDashes - (leftLongestExt.dbLen - leftLongestExt.alnLength + 1), k = 0; k < leftLongestExt.dbLen; i++, k++) {
                        targetAligned[i] = leftLongestExtSeq[k];
                    }      


                    if (outputFileOver.is_open()) {
                        outputFileOver << targetAligned << "\n";
                    }
                }
            }
#endif

        // Iterate through all candidates and compare to the longest
        for (unsigned int leftRes = 0; leftRes < mgeCandiLeft.size(); leftRes++){
            // now retrieve the other candidate extensions and get their sequences

            Matcher::result_t lefty = mgeCandiLeft[leftRes];

            unsigned int leftExtCandiId = sequenceDbr->getId(lefty.dbKey);
            unsigned int leftExtCandiLen = sequenceDbr->getSeqLen(leftExtCandiId);

            if ( leftLongestExtId == leftExtCandiId ){
                continue;
            }

            char *leftExtCandiSequ = sequenceDbr->getData(leftExtCandiId, thread_idx);
            std::string leftExtCandiSeq;

            //bool isBadCandi = false;
            //std::vector<long double> expSim;
            if (lefty.isRevToAlignment) {
                // Convert the reversed fragment to std::string
                char *leftExtCandiSeqTmp = getNuclRevFragment(leftExtCandiSequ, leftExtCandiLen, (NucleotideMatrix *)subMat);
                leftExtCandiSeq = std::string(leftExtCandiSeqTmp, leftExtCandiLen);
                delete[] leftExtCandiSeqTmp;
                // Here we calculate the epected similarity of the extension
                //expSim = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNucRev, seqErrMatch, false);
            } else {
                // Directly convert the raw sequence to std::string
                leftExtCandiSeq = std::string(leftExtCandiSequ, leftExtCandiLen);
                // Here we calculate the epected similarity of the extension
                //expSim = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNuc, seqErrMatch, false);
            }

            // calculate the sequence identity
            unsigned int othersExtLen = lefty.dbLen - lefty.alnLength;
            int idCnt = 0;
            int idRyCnt = 0;
            for (unsigned int leftPos = 0; leftPos < othersExtLen; ++leftPos) {
                idCnt += (leftLongestExtSeq[leftLongestExt.dbStartPos - 1 - leftPos] == leftExtCandiSeq[lefty.dbStartPos - 1 - leftPos]) ? 1 : 0;
                idRyCnt += (ryMap[leftLongestExtSeq[leftLongestExt.dbStartPos - 1 - leftPos]] == ryMap[leftExtCandiSeq[lefty.dbStartPos - 1 - leftPos]]) ? 1 : 0;
            }

            //float seqId = static_cast<float>(idCnt) / othersExtLen;
            float rySeqId = static_cast<float>(idRyCnt) / othersExtLen;
            
            // check for number of C-T and G-A mismatches:
            unsigned int ryMM = idRyCnt - idCnt;
            unsigned int ryMMThr = othersExtLen * 0.1; //allowing for 10% mismatch rate

            // //if ( rySeqId < rymerThresh ){
            // //if ( rySeqId < rymerThresh && ( (othersExtLen <= 10 && ryMM <= 2) || ryMM <= ryMMThr )){
            // if ( rySeqId < rymerThresh || ( (othersExtLen < 10 && ryMM <= 2) || seqId < 0.95 )){
            //     mgeFoundLeft = true;
            //     break;
            // } 

            // if (rySeqId < rymerThresh) {
            //     mgeFoundLeft = true;
            //     break;
            // } else if (othersExtLen <= 10 && ryMM > 1) {
            //     mgeFoundLeft = true;
            //     break;
            // }

            if (othersExtLen <= 10 && ryMM > 1) {
                mgeFoundLeft = true;
                break;
            } else if (rySeqId < rymerThresh) {
                mgeFoundLeft = true;
                break;
            } else if ( ryMM > ryMMThr ){
                mgeFoundLeft = true;
                break;  
            }

#ifdef DEBUGCORR
        if ( queryKey == 4631784 ){
        //if ( true ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '-');

                for (int i = leadingDashes - (lefty.dbLen - lefty.alnLength + 1), k = 0; k < lefty.dbLen; i++, k++) {
                    targetAligned[i] = leftExtCandiSeq[k];
                }                       

                // Append the loop output to outputLine
                targetAligned += "\t";
                targetAligned += std::to_string(lefty.seqId);
                targetAligned += "\t";
                // targetAligned += std::to_string(insideSeqId);
                // targetAligned += "\t";
                targetAligned += std::to_string(lefty.rySeqId);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[0]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[1]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[2]);
                targetAligned += "\n";

                if (outputFileOver.is_open()) {
                    outputFileOver << targetAligned;
                }
            }
        }
#endif
        }
    }

    return std::make_pair(mgeFoundLeft,mgeFoundRight);
    // MGE marker end
}













std::pair<bool, bool> mgeFinder(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat, const std::vector<diNucleotideProb> & subDeamDiNuc, const std::vector<diNucleotideProb> & subDeamDiNucRev, const diNucleotideProb & seqErrMatch){

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

// MGE marker
    bool mgeFoundRight = false;
    bool mgeFoundLeft = false;

    std::vector<Matcher::result_t> mgeCandiRight;
    std::vector<Matcher::result_t> mgeCandiLeft;

    // define thresholds to include extension candidates for the mge-identifier
    // float alnSeqIdThr = par.seqIdThr;
    //float alnSeqIdThr = 0.999;
    //unsigned int extLen = 150;

    //set threshold when to declare an mge
    // float mgeSeqId = 0.9;
    // float mgeRySeqId = 0.99;

    //set number of bases in overlap to comapare after the alignment ends
    //unsigned int numBaseCompare = 50;

    for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

        Matcher::result_t res = alignments[alnIdx];
        size_t dbKey = res.dbKey;
        unsigned int resId = sequenceDbr->getId(res.dbKey);
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);


        // distinguish between right and left here

        unsigned int dbStartPos = res.dbStartPos;
        unsigned int dbEndPos = res.dbEndPos;
        unsigned int qStartPos = res.qStartPos;
        unsigned int qEndPos = res.qEndPos;
        unsigned int targetId = sequenceDbr->getId(res.dbKey);
        if (targetId == UINT_MAX) {
            Debug(Debug::ERROR) << "Could not find targetId  " << res.dbKey
                                << " in database " << sequenceDbr->getDataFileName() << "\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId);

        // TODO: Add min overlap len and max target seq len
        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){

            char *resSeq = sequenceDbr->getData(resId, thread_idx); 
            std::string resSeqStr;
            float ryId;

            float rymerThresh = par.correctionThresholdRySeqId;
            //float rymerThresh = 0.95;
            if ( res.alnLength <= 100){
                rymerThresh = (static_cast<float>(res.alnLength) - 1) / static_cast<float>(res.alnLength);
                rymerThresh = std::floor(rymerThresh * 1000) / 1000;
            }

            if ( res.isRevToAlignment ) {
                char *resSeqTmp = getNuclRevFragment(resSeq, res.dbLen, (NucleotideMatrix *) subMat);
                ryId = getRYSeqId(res, querySeq, resSeqTmp, ryMap);
                res.rySeqId = ryId;
                delete[] resSeqTmp;
            }
            else {
                ryId = getRYSeqId(res, querySeq, resSeq, ryMap);
                res.rySeqId = ryId;
            }

            if ( ryId >= rymerThresh ){
                // count only the left and right extensions which extensions should be compared
                //if (res.seqId > alnSeqIdThr && (res.dbLen - res.alnLength <= extLen)){
                //if (res.seqId >= alnSeqIdThr){
                if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                    // right extension
                    mgeCandiRight.push_back(res);
                }
                else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                    // left extension
                    mgeCandiLeft.push_back(res);
                }
                //}
            }
        }
    }


#ifdef DEBUGCORR

        int leadingDashes = 30;
        if ( queryKey == 4631784 ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string headers = "sequences\toverlap_similarity\toverlap_rySimilarity\tsimilarity_to_longest_extension\trySimilarity_to_longest_extension\texpected_similarity";

                std::string qeueryDebug;
                // Construct qeueryDebug from querySeq
                for (int i = 0; i < querySeqLen; i++) {
                    qeueryDebug += querySeq[i];
                }

                // Create a string of 100 "-" characters
                std::string dashes(leadingDashes, '*');

                // Add dashes before and after qeueryDebug
                qeueryDebug = dashes + qeueryDebug + dashes;

                if (outputFileOver.is_open()) {
                    outputFileOver << headers << "\n";
                    outputFileOver << qeueryDebug << "\n";
                }
            }
        }
#endif

    // Now we iterate through the left and right extensions
    // Get the longest extension and align all remaining candidates
    // right extensions
    if (mgeCandiRight.size() > 1){
        
        // Get extension candidate with the longest extension 
        Matcher::result_t rightLongestExt = mgeCandiRight.front();
        for (unsigned int rightCandis = 0; rightCandis < mgeCandiRight.size(); rightCandis++){
            rightLongestExt = (mgeCandiRight[rightCandis].dbLen - mgeCandiRight[rightCandis].alnLength >= rightLongestExt.dbLen - rightLongestExt.alnLength) ? mgeCandiRight[rightCandis] : rightLongestExt;
        }

        // now retrieve the sequence of the longest extension in the correct orientation 
        unsigned int rightLongestExtId = sequenceDbr->getId(rightLongestExt.dbKey);
        unsigned int rightLongestExtLen = sequenceDbr->getSeqLen(rightLongestExtId);

        char *rightLongestExtSequ = sequenceDbr->getData(rightLongestExtId, thread_idx);
        std::string rightLongestExtSeq;

        if ( rightLongestExt.isRevToAlignment ) {
            char *rightLongestExtSeqTmp = getNuclRevFragment(rightLongestExtSequ, rightLongestExtLen, (NucleotideMatrix *) subMat);
            rightLongestExtSeq = std::string(rightLongestExtSeqTmp, rightLongestExtLen);
            delete[] rightLongestExtSeqTmp;
        }
        else {
            rightLongestExtSeq = std::string(rightLongestExtSequ, rightLongestExtLen);
        }
        
#ifdef DEBUGCORR
            if ( queryKey == 4631784 ){
            // DEBUG CORRECTION
                #pragma omp critical
                {
                    std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '+');

                    for (int i = leadingDashes + (querySeqLen - rightLongestExt.alnLength), k = 0; k < rightLongestExt.dbLen; i++, k++) {
                            targetAligned[i] = rightLongestExtSeq[k];
                    }      

                    if (outputFileOver.is_open()) {
                        outputFileOver << targetAligned << "\n";
                    }
                }
            }
#endif

        // Iterate through all candidates and compare to the longest
        for (unsigned int rightRes = 0; rightRes < mgeCandiRight.size(); rightRes++){
            // now retrieve the other candidate extensions and get their sequences

            Matcher::result_t righty = mgeCandiRight[rightRes];

            unsigned int rightExtCandiId = sequenceDbr->getId(righty.dbKey);
            unsigned int rightExtCandiLen = sequenceDbr->getSeqLen(rightExtCandiId);

            if ( rightExtCandiId == rightLongestExtId ){
                continue;
            }

            char *rightExtCandiSequ = sequenceDbr->getData(rightExtCandiId, thread_idx);
            std::string rightExtCandiSeq;

            bool isBadCandi = false;
            std::vector<long double> expSim;
            if (righty.isRevToAlignment) {
                // Convert the reversed fragment to std::string
                char *rightExtCandiSeqTmp = getNuclRevFragment(rightExtCandiSequ, rightExtCandiLen, (NucleotideMatrix *)subMat);
                rightExtCandiSeq = std::string(rightExtCandiSeqTmp, rightExtCandiLen);
                delete[] rightExtCandiSeqTmp;
                // Here we calculate the epected similarity of the extension
                expSim = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNucRev, seqErrMatch, true);
            } else {
                // Directly convert the raw sequence to std::string
                rightExtCandiSeq = std::string(rightExtCandiSequ, rightExtCandiLen);
                // Here we calculate the epected similarity of the extension
                expSim = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNuc, seqErrMatch, true);
            }

#ifdef DEBUGCORR
        if ( queryKey == 4631784 ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '-');
            
                for (int i = leadingDashes + (querySeqLen - righty.alnLength), k = 0; k < righty.dbLen; i++, k++) {
                    targetAligned[i] = rightExtCandiSeq[k];
                }

                // Append the loop output to outputLine
                targetAligned += "\t";
                targetAligned += std::to_string(righty.seqId);
                targetAligned += "\t";
                // targetAligned += std::to_string(insideSeqId);
                // targetAligned += "\t";
                targetAligned += std::to_string(righty.rySeqId);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[0]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[1]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[2]);
                targetAligned += "\n";

                if (outputFileOver.is_open()) {
                    outputFileOver << targetAligned;
                }
            }
        }
#endif

            if ( isBadCandi ){
                mgeFoundRight = true;
                break;
            }
                // unsigned int othersExtLen = ( righty.dbLen - righty.alnLength < numBaseCompare ) ? righty.dbLen - righty.alnLength : numBaseCompare;
                // int idCnt = 0;
                // int idRyCnt = 0;
                // for (unsigned int rightPos = 0; rightPos < othersExtLen; ++rightPos) {
                //     idCnt += (rightLongestExtSeq[rightLongestExt.dbEndPos + rightPos] == rightExtCandiSeq[righty.dbEndPos + rightPos]) ? 1 : 0;
                //     idRyCnt += (ryMap[rightLongestExtSeq[rightLongestExt.dbEndPos + rightPos]] == ryMap[rightExtCandiSeq[righty.dbEndPos + rightPos]]) ? 1 : 0;
                // }

                // float seqId = static_cast<float>(idCnt) / othersExtLen;
                // float rySeqId = static_cast<float>(idRyCnt) / othersExtLen;

                // if ( seqId < mgeSeqId || rySeqId < mgeRySeqId ){
                // //if ( seqId < mgeSeqId ){
                //     // We believe that we found an MGE so do not allow to extend to right at all!!
                //     // std::cerr << "mge RIGHT found:\t" << std::endl;
                //     // std::cerr << "main query\t" << querySeq;
                //     // std::cerr << "rightLongestExtSeq\t" << rightLongestExtSeq << std::endl;
                //     // std::cerr << "rightExtCandiSeq\t" << rightExtCandiSeq << std::endl;
                //     // std::cerr << "seqId:\t" << seqId << std::endl;
                //     // std::cerr << "rySeqId:\t" << rySeqId << std::endl;
                //     // std::cerr << "Compare len:\t" << othersExtLen << std::endl;
                //     // std::cerr << "rightLongestExt.dbStartPos\t" << rightLongestExt.dbStartPos << std::endl;
                //     mgeFoundRight = true;
                //     break;
                // }
        }
    }

    // left extensions
    if (mgeCandiLeft.size() > 1){
        Matcher::result_t leftLongestExt = mgeCandiLeft.front();

        for (unsigned int leftCandis = 0; leftCandis < mgeCandiLeft.size(); leftCandis++){
            leftLongestExt = (mgeCandiLeft[leftCandis].dbLen - mgeCandiLeft[leftCandis].alnLength >= leftLongestExt.dbLen - leftLongestExt.alnLength) ? mgeCandiLeft[leftCandis] : leftLongestExt;
        }

        // now retrieve the sequence of the longest extension in the correct orientation 
        unsigned int leftLongestExtId = sequenceDbr->getId(leftLongestExt.dbKey);
        unsigned int leftLongestExtLen = sequenceDbr->getSeqLen(leftLongestExtId);

        char *leftLongestExtSequ = sequenceDbr->getData(leftLongestExtId, thread_idx);
        std::string leftLongestExtSeq;

        if ( leftLongestExt.isRevToAlignment ) {
            char *leftLongestExtSeqTmp = getNuclRevFragment(leftLongestExtSequ, leftLongestExtLen, (NucleotideMatrix *) subMat);
            leftLongestExtSeq = std::string(leftLongestExtSeqTmp, leftLongestExtLen);
            delete[] leftLongestExtSeqTmp;
        }
        else {
            leftLongestExtSeq = std::string(leftLongestExtSequ, leftLongestExtLen);
        }

#ifdef DEBUGCORR
            if ( queryKey == 4631784 ){
            // DEBUG CORRECTION
                #pragma omp critical
                {
                    std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '+');

                    for (int i = leadingDashes - (leftLongestExt.dbLen - leftLongestExt.alnLength + 1), k = 0; k < leftLongestExt.dbLen; i++, k++) {
                        targetAligned[i] = leftLongestExtSeq[k];
                    }      


                    if (outputFileOver.is_open()) {
                        outputFileOver << targetAligned << "\n";
                    }
                }
            }
#endif

        // Iterate through all candidates and compare to the longest
        for (unsigned int leftRes = 0; leftRes < mgeCandiLeft.size(); leftRes++){
            // now retrieve the other candidate extensions and get their sequences

            Matcher::result_t lefty = mgeCandiLeft[leftRes];

            unsigned int leftExtCandiId = sequenceDbr->getId(lefty.dbKey);
            unsigned int leftExtCandiLen = sequenceDbr->getSeqLen(leftExtCandiId);

            if ( leftLongestExtId == leftExtCandiId ){
                continue;
            }

            char *leftExtCandiSequ = sequenceDbr->getData(leftExtCandiId, thread_idx);
            std::string leftExtCandiSeq;

            bool isBadCandi = false;
            std::vector<long double> expSim;
            if (lefty.isRevToAlignment) {
                // Convert the reversed fragment to std::string
                char *leftExtCandiSeqTmp = getNuclRevFragment(leftExtCandiSequ, leftExtCandiLen, (NucleotideMatrix *)subMat);
                leftExtCandiSeq = std::string(leftExtCandiSeqTmp, leftExtCandiLen);
                delete[] leftExtCandiSeqTmp;
                // Here we calculate the epected similarity of the extension
                expSim = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNucRev, seqErrMatch, false);
            } else {
                // Directly convert the raw sequence to std::string
                leftExtCandiSeq = std::string(leftExtCandiSequ, leftExtCandiLen);
                // Here we calculate the epected similarity of the extension
                expSim = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNuc, seqErrMatch, false);
            }

#ifdef DEBUGCORR
        if ( queryKey == 4631784 ){
        //if ( true ){
        // DEBUG CORRECTION
            #pragma omp critical
            {
                std::string targetAligned(querySeqLen+leadingDashes+leadingDashes, '-');

                for (int i = leadingDashes - (lefty.dbLen - lefty.alnLength + 1), k = 0; k < lefty.dbLen; i++, k++) {
                    targetAligned[i] = leftExtCandiSeq[k];
                }                       

                // Append the loop output to outputLine
                targetAligned += "\t";
                targetAligned += std::to_string(lefty.seqId);
                targetAligned += "\t";
                // targetAligned += std::to_string(insideSeqId);
                // targetAligned += "\t";
                targetAligned += std::to_string(lefty.rySeqId);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[0]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[1]);
                targetAligned += "\t";
                targetAligned += std::to_string(expSim[2]);
                targetAligned += "\n";

                if (outputFileOver.is_open()) {
                    outputFileOver << targetAligned;
                }
            }
        }
#endif

            if ( isBadCandi ){
                mgeFoundLeft = true;
                break;
            }
            // // calculate the sequence identity
            // unsigned int othersExtLen = ( lefty.dbLen - lefty.alnLength < numBaseCompare ) ? lefty.dbLen - lefty.alnLength : numBaseCompare;
            // int idCnt = 0;
            // int idRyCnt = 0;
            // for (unsigned int leftPos = 0; leftPos < othersExtLen; ++leftPos) {
            //     // std::cerr << "position\t" << leftPos << std::endl;
            //     // std::cerr << "leftLongestExt.dbStartPos\t" << leftLongestExt.dbStartPos << std::endl;
            //     // std::cerr << "leftLongestExt.dbEndPos\t" << leftLongestExt.dbEndPos << std::endl;
            //     // std::cerr << "lefty.dbStartPos\t" << lefty.dbStartPos << std::endl;
            //     // std::cerr << "lefty.dbEndPos\t" << lefty.dbEndPos << std::endl;
            //     // std::cerr << "leftLongestExtSeq\t" << leftLongestExtSeq << std::endl;
            //     // std::cerr << "leftExtCandiSeq\t" << leftExtCandiSeq << std::endl;
            //     idCnt += (leftLongestExtSeq[leftLongestExt.dbStartPos - 1 - leftPos] == leftExtCandiSeq[lefty.dbStartPos - 1 - leftPos]) ? 1 : 0;
            //     idRyCnt += (ryMap[leftLongestExtSeq[leftLongestExt.dbStartPos - 1 - leftPos]] == ryMap[leftExtCandiSeq[lefty.dbStartPos - 1 - leftPos]]) ? 1 : 0;
            // }

            // float seqId = static_cast<float>(idCnt) / othersExtLen;
            // float rySeqId = static_cast<float>(idRyCnt) / othersExtLen;

            // if ( seqId < mgeSeqId || rySeqId < mgeRySeqId ){
            // //if ( seqId < mgeSeqId ){
            //     // We believe that we found an MGE so do not allow to extend to LEFT at all!!
            //     // std::cerr << "mge LEFT found:\t" << std::endl;
            //     // std::cerr << "main query\t" << querySeq;
            //     // std::cerr << "leftLongestExtSeq\t" << leftLongestExtSeq << std::endl;
            //     // std::cerr << "leftExtCandiSeq\t" << leftExtCandiSeq << std::endl;
            //     // std::cerr << "seqId:\t" << seqId << std::endl;
            //     // std::cerr << "rySeqId:\t" << rySeqId << std::endl;
            //     // std::cerr << "Compare len:\t" << othersExtLen << std::endl;
            //     // std::cerr << "leftLongestExt.dbStartPos\t" << leftLongestExt.dbStartPos << std::endl;
            //     // std::cerr << std::endl; 
            //     mgeFoundLeft = true;
            //     break;
            // }
        }
    }

    return std::make_pair(mgeFoundLeft,mgeFoundRight);
    // MGE marker end

}




























std::vector<long double> calcLikelihoodCorrection(const Matcher::result_t & rightLongestCandi, const Matcher::result_t & candidate, const std::string & querySeq, const std::string & targetSeq, const std::vector<diNucleotideProb> & subDeamDiNuc, const diNucleotideProb & seqErrMatch, bool extSide)
{
    std::vector<float> baseFreqs = { 0.23554, 0.26446, 0.26446, 0.23554 };

    std::unordered_map<char, int> ryMap = {
    {'A', 0},
    {'C', 1},
    {'G', 0},
    {'T', 1}};

    std::unordered_map<char, int> nucleotideMap = {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3}};

    unsigned int countMatch = 0;
    unsigned int countMismatch = 0;

    unsigned int countRyMatch = 0;
    unsigned int countRyMismatch = 0;

    // Here we start calculating the likelihood for a >>> S I N G L E <<< alignment
    long double likMod = -std::numeric_limits<long double>::infinity(); // log(0)

    std::vector<diNucleotideProb> subdeam_lookup(candidate.dbLen);

    // Create lookup vector
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[i];
    }
    for (size_t i = 5; i < candidate.dbLen - 5; ++i) {
        subdeam_lookup[i] = subDeamDiNuc[5];
    }
    for (size_t i = 0; i < 5; ++i) {
        subdeam_lookup[candidate.dbLen - 5 + i] = subDeamDiNuc[ 6 + i];
    }


// Actually I think you don't even need to have the sequences available here...
// Because this calculation is GIVEN THAT THERE IS A RYMER MATCH
// However what we need is the deamination pattern at the particular positions
    unsigned int toffset = 0;
    unsigned int qoffset = 0;
    const std::string rymerBase = "1010";

    //unsigned int longest = (rightLongestCandi.dbLen > candidate.dbLen) ? rightLongestCandi.dbLen : candidate.dbLen;

    unsigned int minExt = candidate.dbLen - candidate.alnLength;
    unsigned int minALn = (rightLongestCandi.alnLength < candidate.alnLength) ? rightLongestCandi.alnLength : candidate.alnLength;
    unsigned int minTotal = minExt + minALn;


    // right extension; actual bases only needed for the observed sequence identity
    if ( extSide == true ){

        if (candidate.qStartPos < rightLongestCandi.qStartPos) {
            toffset = rightLongestCandi.qStartPos - candidate.qStartPos;
        } else if (rightLongestCandi.qStartPos < candidate.qStartPos) {
            qoffset = candidate.qStartPos - rightLongestCandi.qStartPos;
        }

        for (unsigned int t = 0; t < minTotal; t++) {
            
            unsigned int qIndex = rightLongestCandi.dbStartPos + qoffset + t;
            unsigned int tIndex = candidate.dbStartPos + toffset + t;

            int qBase = nucleotideMap[querySeq[qIndex]];
            int tBase = nucleotideMap[targetSeq[tIndex]];
            int qBaseRy = ryMap[querySeq[qIndex]];
            int tBaseRy = ryMap[targetSeq[tIndex]];

            countMatch += (tBase == qBase) ? 1 : 0;
            countMismatch += (tBase != qBase) ? 1 : 0;
            countRyMatch += (tBaseRy == qBaseRy) ? 1 : 0;
            countRyMismatch += (tBaseRy != qBaseRy) ? 1 : 0;
        
            long double mmLik = -std::numeric_limits<long double>::infinity(); // log(0)
            diNucleotideProb tProbs = subdeam_lookup[candidate.dbStartPos + toffset + t];

            for (int qIt = 0; qIt < 4; qIt++){
                long double likeli = 0;
                for (int tIt = 0; tIt < 4; tIt++){
                    if ( rymerBase[qIt] != rymerBase[tIt] ){
                        continue;
                    }
                    long double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qIt][tIt]);
                    likeli += (std::log(baseFreqs[qIt]) + std::log(match_lik));  
                }
                mmLik = oplusnatl(mmLik, likeli);
            }
            likMod = oplusnatl(likMod, mmLik);
        }
    }
    // left extension; actual bases only needed for the observed sequence identity
    else {
        // Left extension
        if (candidate.qEndPos > rightLongestCandi.qEndPos) {
            toffset = candidate.qEndPos - rightLongestCandi.qEndPos;
        } else if (rightLongestCandi.qEndPos > candidate.qEndPos) {
            qoffset = rightLongestCandi.qEndPos - candidate.qEndPos;
        }

        for (unsigned int t = 0; t < minTotal; t++) {
            unsigned int qIndex = rightLongestCandi.dbEndPos - t - qoffset;
            unsigned int tIndex = candidate.dbEndPos - t - toffset;

            //if (qIndex < querySeq.size() && tIndex < targetSeq.size()) {
                int qBase = nucleotideMap[querySeq[qIndex]];
                int tBase = nucleotideMap[targetSeq[tIndex]];
                int qBaseRy = ryMap[querySeq[qIndex]];
                int tBaseRy = ryMap[targetSeq[tIndex]];

                countMatch += (tBase == qBase) ? 1 : 0;
                countMismatch += (tBase != qBase) ? 1 : 0;
                countRyMatch += (tBaseRy == qBaseRy) ? 1 : 0;
                countRyMismatch += (tBaseRy != qBaseRy) ? 1 : 0;
        
                //long double mmLik = 0;
                diNucleotideProb tProbs = subdeam_lookup[candidate.dbEndPos - t - toffset];

                // Initialize log space likelihood
                long double mmLik = -std::numeric_limits<long double>::infinity(); // log(0)

                // Iterate over possible states
                for (int qIt = 0; qIt < 4; qIt++) {
                    long double likeli = 0;
                    for (int tIt = 0; tIt < 4; tIt++) {
                        if (rymerBase[qIt] != rymerBase[tIt]) {
                            continue;
                        }
                        long double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qIt][tIt]);
                        likeli += std::log(baseFreqs[qIt]) + std::log(match_lik);  
                    }
                    mmLik = oplusnatl(mmLik, likeli);
                }
                likMod = oplusnatl(likMod, mmLik);
            }
    }

    unsigned int overlapLen = countMatch + countMismatch;
    // Adjust for the length
    long double adjustment = std::log(overlapLen);
    likMod = likMod - adjustment;

    // Convert back from log space
    likMod = 1 - expl(likMod);

    // Round the expectation value down
    // long double likRounded = 0;
    // for (unsigned int k = 0; k < overlapLen; ++k) {
    //     long double lowerBound = static_cast<float>(overlapLen - k - 1) / overlapLen;
    //     if (likMod >= lowerBound) {
    //         likRounded = lowerBound;
    //         break;
    //     }
    // }

    // get sequence identity 
    long double obsSeqId = static_cast<double>(countMatch)/static_cast<double>(overlapLen);
    long double obsRySeqId = static_cast<double>(countRyMatch)/static_cast<double>(overlapLen);

    std::vector<long double> result = {obsSeqId, obsRySeqId, likMod};

    return result;



    // this is only for the right extension for now
    // Now find the position where rightLongestCandi and candidate are overlapping first
    // unsigned int toffset = 0;
    // unsigned int qoffset = 0;
    // bool candHasLonger;
    // unsigned int startComp = 0;
    // if (candidate.qStartPos < rightLongestCandi.qStartPos){
    //     toffset = rightLongestCandi.qStartPos - candidate.qStartPos;
    //     startComp = toffset;
    // }
    // else if (candidate.qStartPos > rightLongestCandi.qStartPos){
    //     qoffset = candidate.qStartPos - rightLongestCandi.qStartPos;
    // }

    // //We iterate over the full length of the longest extension

    // for ( unsigned int t = toffset; t < candidate.dbLen; t++)
    // {
    //     long double mmLik = 0;

    //     diNucleotideProb tProbs;
    //     tProbs = subdeam_lookup[candidate.qStartPos + t];

    //     // first sequence is reference
    //     int qBase = nucleotideMap[querySeq[rightLongestCandi.dbStartPos + qoffset + t - toffset]];
    //     int tBase = nucleotideMap[targetSeq[candidate.dbStartPos + qoffset + t]];
    //     int qBaseRy = ryMap[querySeq[rightLongestCandi.dbStartPos + qoffset + t - toffset]];
    //     int tBaseRy = ryMap[targetSeq[candidate.dbStartPos + qoffset + t]];
    // }
    // // Iterate over the full length of the longest extension
    // for (unsigned int t = offset; t < candidate.dbLen; t++) {
    //     long double mmLik = 0;

    //     diNucleotideProb tProbs = subdeam_lookup[candidate.qStartPos + t];

    //     // First sequence is reference
    //     int qBase = nucleotideMap[querySeq[rightLongestCandi.dbStartPos + t - offset]];
    //     int tBase = nucleotideMap[targetSeq[candidate.dbStartPos + t]];
    //     int qBaseRy = ryMap[querySeq[rightLongestCandi.dbStartPos + t - offset]];
    //     int tBaseRy = ryMap[targetSeq[candidate.dbStartPos + t]];

    //     for (int qIt = 0; qIt < 4; qIt++){

    //         long double likeli = 0;
    //         for (int tIt = 0; tIt < 4; tIt++){

    //             if ( rymerBase[qBase] != rymerBase[tIt] ){
    //                 continue;
    //             }

    //             double match_lik = std::max(static_cast<long double>(SMOOTHING_VALUE), tProbs.p[qIt][tIt]);
    //             likeli += (std::log(baseFreqs[qIt]) + std::log(match_lik));  
    //         }
    //         mmLik += std::exp(likeli);
    //     }
    // }


    //if ( candidate.seqId < likFin )
    //if ( qKey == 4631784 )
    //if ( obsRySeqId < 0.999 && obsSeqId < likRounded )
    if ( false )
    //if ( candidate.dbKey == 66487)
    //if ( randAln > likMod )
    {
        //std::cerr << "\n";
        //std::cerr << "q_full\t" << querySeq << std::endl;
        //std::cerr << "t_full\t" << targetSeq << std::endl;
        //std::cerr << "q_overlap\t" << queryOverlap << std::endl;
        //std::cerr << "t_overlap\t" << targetOverlap << std::endl;
        //std::cerr << "Aligned seqId:\t" << candidate.seqId << std::endl;
        //std::cerr << "Observed seqId:\t" << obsSeqId << std::endl;
        //std::cerr << "Extension orientation:\t" << extSide << std::endl;
        //std::cerr << "subSeqId:\t" << insideSeqId << std::endl;
        //std::cerr << "rySeqId:\t" << candidate.rySeqId << std::endl;
        //std::cerr << "likMod:\t" << likAlnLen << std::endl;
        //std::cerr << "likMod1:\t" << likMod1 << std::endl;
        //std::cerr << "likRounded:\t" << likRounded << std::endl;
        //std::cerr << "likMod2:\t" << likMod2 << std::endl;
        //std::cerr << "minExt:\t" << minExt << std::endl;
        //std::cerr << "sumMatch:\t" << sumMatches << std::endl;
        //std::cerr << "likFin:\t" << likFin << std::endl;
        //std::cerr << "M and MM:\t" << countMatch << "\t" << countMismatch << std::endl;
        //std::cerr << "TdbKey,Tlen:\t" << candidate.dbKey << "\t" << candidate.dbLen << std::endl;
        //std::cerr << "query length:\t" << candidate.qLen << std::endl;
        // std::cerr << "ratio space" << std::endl;
        // std::cerr << "ratio lik/(lik+random):\t" << ratioLog << std::endl;
        //std::cerr << "alnLen:\t" << candidate.alnLength << std::endl;
        //std::cerr << "\n";
    }    

    // // sequences diverge to much
    // if ( obsSeqId < likRounded || obsRySeqId <= rySeqIdThr ){
    //     return true;
    // }
    // else{ // sequences are within expected bonds
    //     return false;
    // }
    // return false;

} // END loop through each base


// DEBUG TOOLS TO DELETE
// void printSubstitutionRatesVector(const std::vector<substitutionRates>& sub3p) {
//     for (const auto& rates : sub3p) {
//         for (int i = 0; i < 12; ++i) {
//             std::cerr << rates.s[i] << " ";
//         }
//         std::cerr << std::endl; // Print a new line after each substitutionRates struct
//     }
// }

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





