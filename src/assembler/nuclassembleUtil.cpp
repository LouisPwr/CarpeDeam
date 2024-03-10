#include "nuclassembleUtil.h"

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

void getSeqErrorProf(diNucleotideProb & seqErrMatch, diNucleotideProb & seqErrMis, long double err)
{
    for (int orig = 0; orig < 4; orig++)
    {
        for (int obs = 0; obs < 4; obs++)
        {
            seqErrMatch.p[orig][obs] = (orig == obs) ? 1 - err : err/3;
        }
    }

    for (int orig = 0; orig < 4; orig++)
    {
        for (int obs = 0; obs < 4; obs++)
        {
            seqErrMis.p[orig][obs] = (( orig == 1 && obs == 3 ) || ( orig == 2 && obs == 0 )) ? 1 - err : err/3;
        }
    }
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


// calculate number of matches/mismatches for adapted bayesian extension in nuclassembleresult2.cpp
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
        unsigned int subDeamLen = subDeamDiNuc.size();
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);

        // // right overlap
        // if ( rightStart ){
        //     //std::cerr << "right start" << std::endl;
        //     if ( pos < 5 ){
        //         matchLik = subDeamDiNuc[pos].p[qBase][tBase];
        //     }
        //     else if ( pos >= res.dbLen - 5 ){
        //         matchLik = subDeamDiNuc[subDeamLen - (res.alnLength - pos)].p[qBase][tBase];
        //     }
        //     else{
        //         matchLik = subDeamDiNuc[5].p[qBase][tBase];
        //     }
        // }
        // //left overlap
        // else if ( leftStart ){
        //     //std::cerr << "left start" << std::endl;
        //     // Iterating over the alignment
        //     if ( res.dbStartPos + pos < 5 ){
        //         matchLik = subDeamDiNuc[res.dbStartPos + pos].p[qBase][tBase];
        //     }
        //     else if ( pos >= res.alnLength - 5 ){
        //         matchLik = subDeamDiNuc[subDeamLen - (res.alnLength - pos)].p[qBase][tBase];
        //     }
        //     else{
        //         matchLik = subDeamDiNuc[5].p[qBase][tBase];
        //     }
        // }
        matchLik = subDeamDiNuc[5].p[qBase][tBase];

        if ( dimerCT && matchLik > 0 ){
            mCT += deamMatches(res, scoreAln, matchLik);
            // if ( matchLik > 0.01 ){
            // float matches = ((static_cast<float>(scoreAln) + 3.0f*res.alnLength) / 5.0f) + mCT + mGA; // TODO: use match score from nucleotide.out
            // float matchesWO = (static_cast<float>(scoreAln) + 3.0f*res.alnLength) / 5.0f; // TODO: use match score from nucleotide.out matrix
            // float mOld = res.seqId * res.alnLength ;
            // std::cerr << "matches new " << matches << std::endl;
            // std::cerr << "matches without deam " << matchesWO << std::endl;
            // std::cerr << "matches conservative " << mOld << std::endl; 
            // std::cerr << "aln len " << res.alnLength << std::endl;
            // std::cerr << "score of aln " << scoreAln << std::endl;
            // std::cerr << "matchLik " << matchLik << std::endl;
            // }
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


void calc_likelihood(scorePerRes & scoredRes, char* querySeq, const char* targetSeq, std::vector<diNucleotideProb> & subDeamDiNuc, unsigned int maxAln, bool isRightOverlap, float randAlnPenal, diNucleotideProb & seqErrMatch, diNucleotideProb & seqErrMis)
{
    unsigned int countMatch = 0;
    unsigned int countMismatch = 0;

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
        double likMod = 0.0; //used to save sum of all log_lik_r1

        std::unordered_map<char, int> nucleotideMap = {
        {'A', 0},
        {'C', 1},
        {'G', 2},
        {'T', 3}};


    //You either have a match or a mismatch
    for ( unsigned int t = 0; t < scoredRes.r.alnLength; t++ )
    {   
        double lik = 0;
        unsigned int subDeamLen = subDeamDiNuc.size();

        diNucleotideProb tProbs;

        // for both left and right extension as the matrices were rotated accordingly before this function is called

        // We know that the overlap is always at least k-long (i.e. 22)
        // This can be changed to a better/more efficient solution once we know it works!
        if ( isRightOverlap )
        {
            if ( t < 5 )
            {
                tProbs = subDeamDiNuc[t];
            }
            else if ( t >= scoredRes.r.dbLen - 5 )
            {
                tProbs = subDeamDiNuc[subDeamLen - (scoredRes.r.alnLength - t)];
            }
            else
            {
                tProbs = subDeamDiNuc[5];
            }
        }
        else
        {
            // left Extension
            if ( scoredRes.r.dbStartPos + t < 5 )
            {
                tProbs = subDeamDiNuc[scoredRes.r.dbStartPos + t];
            }
            else if ( t >= scoredRes.r.alnLength - 5 )
            {
                tProbs = subDeamDiNuc[subDeamLen - (scoredRes.r.alnLength - t)];
            }
            else{
                tProbs = subDeamDiNuc[5];
            }
        }

        // MATCH
        if ( querySeq[scoredRes.r.qStartPos + t] == targetSeq[scoredRes.r.dbStartPos + t] )
        {
            for (int query = 0; query<4; query++)
            {   
                int t_base = nucleotideMap[targetSeq[scoredRes.r.dbStartPos + t]];
                double match_lik = tProbs.p[query][t_base];

                lik += match_lik * seqErrMatch.p[query][t_base];

            }
            countMatch += 1;
        }
        // Mismatch
        else
        {
            for (int query = 0; query<4; query++)
            {   
                int t_base = nucleotideMap[targetSeq[scoredRes.r.dbStartPos + t]];
                double mm_lik = 0;

                mm_lik = tProbs.p[query][t_base];
                lik += mm_lik * seqErrMis.p[query][t_base];

            }
            countMismatch += 1;
        }

    likMod += log(lik);

    } // else statment ending of the S or rate bases case senario


    // likelihood w/o penalty
    double likAlnLen = likMod;

    // Penalizing short alignments
    unsigned int excess = maxAln - scoredRes.r.alnLength;
    likMod += ( excess * log(1.0/4.0) );

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

// void createNoDamageMatrix(const std::string& filename) {
//     std::ofstream file(filename);

//     if (!file.is_open()) {
//         std::cerr << "Error opening file: " << filename << std::endl;
//         return;
//     }

//     // Header row
//     file << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";

//     // Matrix rows
//     for (int i = 0; i < 5; ++i) {
//         for (int j = 0; j < 11; ++j) {
//             file << "0.0\t";
//         }
//         file << "0.0\n"; // Last element in the row without a trailing tab
//     }

//     file.close();
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
        const auto& origStruct = sub3p.back(); // Directly access the last element
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

    for (auto rit = sub3p.rbegin(); rit != sub3p.rend(); ++rit) {
        // Initialize new struct
        diNucleotideProb singleDeam;

        // Access the original struct like so
        auto& origStruct = *rit;

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





