#include "nuclassembleUtil.h"
const double SMOOTHING_VALUE = 0.001;

//#define DEBUG_CORR

int mostLikeliBaseRead(const int baseInQuery, const int qIter, const std::vector<countDeamCov> & deamVec, const std::vector<countDeamCov> & countRevs, const std::vector<diNucleotideProb> & subDeamDiNuc, const std::vector<diNucleotideProb> & subDeamDiNucRev, const diNucleotideProb & seqErrMatch)
{
    std::vector<double> baseLikelis; 

    // Initial likelihod of qBase (either A,C,G or T)
    //std::vector<float> baseFreqs = { 0.23554, 0.26446, 0.26446, 0.23554 };
    std::vector<float> baseFreqs = { 0.25, 0.25, 0.25, 0.25 };

    for ( int qBase=0; qBase<4; qBase++ )
    {
        double qBaseLik = 0;
        long double qBaseErr = 0;
        qBaseErr = seqErrMatch.p[qBase][baseInQuery];
        for (int tBase=0; tBase < 4; tBase++)
        {
            // The vector deamVec contains the counts of deaminations per position;
            // Iterating through it to make the amount of multiplications we need to do. 
            for ( unsigned int l = 0; l < subDeamDiNuc.size(); l++)
            {
                double deamPattern = subDeamDiNuc[l].p[qBase][tBase];
                double deamPatternRev = subDeamDiNucRev[l].p[qBase][tBase];
                // this can only be true in case of a mismatch of Pur<->Pyr
                deamPattern = std::max(deamPattern, SMOOTHING_VALUE);
                //(deamPattern == 0.0) ? SMOOTHING_VALUE : deamPattern;
                deamPatternRev = std::max(deamPatternRev, SMOOTHING_VALUE);
                //(deamPatternRev == 0.0) ? SMOOTHING_VALUE : deamPatternRev;

                int covDeam = deamVec[qIter].count[tBase][l];
                int numReverse = countRevs[qIter].count[tBase][l];
                //qBaseLik += (covDeam - numReverse)*( log( baseFreqs[qBase] ) + log(qBaseErr) + log(deamPattern) );
                //qBaseLik += numReverse*( log( baseFreqs[qBase] ) + log(qBaseErr) + log(deamPatternRev) );
                qBaseLik += (covDeam - numReverse)*( log( baseFreqs[qBase] ) + log(qBaseErr) + log(deamPattern) );
                qBaseLik += numReverse*( log( baseFreqs[qBase] ) + log(deamPatternRev) );
            }
        }
        baseLikelis.push_back(qBaseLik);
    }

    //for (int i = 0; i<4; i++)
    //{
    //    std::cerr << "Base:\t" << i << "\tLikeli\t" << baseLikelis[i] << std::endl;
    //
/*     std::cerr << "IdMax:\t" << idxMax << std::endl;
    std::cerr << '\n' << std::endl; */
    auto maxLik = std::max_element(baseLikelis.begin(), baseLikelis.end());
    int idxMax = std::distance(baseLikelis.begin(), maxLik);
    return idxMax;
}



int doCorrection(LocalParameters &par) {

    DBReader<unsigned int> *sequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * alnReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, sequenceDbr->getDbtype());
    resultWriter.open();

    int seqType = sequenceDbr->getDbtype();
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    }

    EvalueComputation evaluer(sequenceDbr->getAminoAcidDBSize(), subMat);

    unsigned char * wasExtended = new unsigned char[sequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+sequenceDbr->getSize(), 0);
    Debug::Progress progress(sequenceDbr->getSize());

    std::string userInput = par.ancientDamagePath; // This should come from the user
    
    std::string high5 = userInput + "5p.prof";
    std::string high3 = userInput + "3p.prof";

#ifdef DEBUG_CORR
    auto currentTime = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(currentTime);

    std::stringstream one;
    one << "coverages_" << time << ".tsv";
    std::string filenameOne = one.str();
    std::ofstream outputFileAll;  // Declare the ofstream object globally if possible
    outputFileAll.open(filenameOne, std::ios::app);  // Open once, append mode

    std::stringstream two;
    two << "overlaps_" << time << ".tsv";
    std::string filenameTwo = two.str();
    std::ofstream outputFileOver;  // Declare the ofstream object globally if possible
    outputFileOver.open(filenameTwo, std::ios::app);  // Open once, append mode

#endif

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        std::vector<Matcher::result_t> alignments;
        alignments.reserve(300);
        bool *useReverse = new bool[sequenceDbr->getSize()];
        std::fill(useReverse, useReverse+sequenceDbr->getSize(), false);

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


        // Preparation of damage related variables and data structures
        std::vector<diNucleotideProb> subDeamDiNuc(11);
        //Substitution rates due to deamination
        std::vector<substitutionRates> sub5p;
        std::vector<substitutionRates> sub3p;

        std::vector<diNucleotideProb> subDeamDiNucRev;

        //std::cerr << "Before initDeamProbabilities" << std::endl;
        initDeamProbabilities(high5, high3, sub5p, sub3p, subDeamDiNuc, subDeamDiNucRev);

        diNucleotideProb seqErrMatch;

        long double seqErrCorrection= 0.01;
        getSeqErrorProf(seqErrMatch, seqErrCorrection);


#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {

            progress.updateProgress();

            unsigned int queryKey = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id, thread_idx);
            unsigned int querySeqLen = sequenceDbr->getSeqLen(id);
            //std::string query(querySeq, querySeqLen); // no /n/0

            char *alnData = alnReader->getDataByDBKey(queryKey, thread_idx);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);

            bool qWasExtended = sequenceDbr->getExtData(id);

            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

                int rawScore = static_cast<int>(evaluer.computeRawScoreFromBitScore(alignments[alnIdx].score) + 0.5);
                float scorePerCol = static_cast<float>(rawScore) / static_cast<float>(alignments[alnIdx].alnLength + 0.5);

                alignments[alnIdx].score = static_cast<int>(scorePerCol*100);

                // here I will also change the target sequence so that all of them are in the right orientation towards the query
                // update: No, I wont, because changing the data is not per se done fast as sequences are written in ht database in each iteration
                if (seqType == Parameters::DBTYPE_NUCLEOTIDES) {
                    if (alignments[alnIdx].qStartPos > alignments[alnIdx].qEndPos) {
                        useReverse[sequenceDbr->getId(alignments[alnIdx].dbKey)] = true;
                        std::swap(alignments[alnIdx].qStartPos, alignments[alnIdx].qEndPos);
                        unsigned int dbStartPos = alignments[alnIdx].dbStartPos;
                        alignments[alnIdx].dbStartPos = alignments[alnIdx].dbLen - alignments[alnIdx].dbEndPos - 1;
                        alignments[alnIdx].dbEndPos= alignments[alnIdx].dbLen - dbStartPos - 1;
                        alignments[alnIdx].isRevToAlignment = true;

                    } else {
                        useReverse[sequenceDbr->getId(alignments[alnIdx].dbKey)] = false;
                        alignments[alnIdx].isRevToAlignment = false;
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Need to change the deams as they change for each target !!! //
            
            // All you need is the coverage for each position and ...
            std::vector<std::vector<unsigned int>> queryCov(querySeqLen, std::vector<unsigned int> (4, 0));
            std::vector<unsigned int> totalCov(querySeqLen, 0);

            // vector of length of query; for each position we have a vector of lenght deamination matrix (=11)
            // marking for each target which deamination pattern it will refer to in the query 
            std::vector<countDeamCov> deamVec(querySeqLen);
            std::vector<countDeamCov> revCount(querySeqLen);

            // Initialize the counts with zero for each element in the vector 
            for (unsigned int i = 0; i < querySeqLen; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 11; k++) {
                        deamVec[i].count[j][k] = 0;
                        revCount[i].count[j][k] = 0;
                    }
                }
            }


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

                    float rymerThresh = par.correctionThreshold;
                    //float rymerThresh = 0.95;
                    if ( res.alnLength <= 100){
                        rymerThresh = (static_cast<float>(res.alnLength) - 1) / static_cast<float>(res.alnLength);
                        rymerThresh = std::floor(rymerThresh * 1000) / 1000;
                    }

                    if ( res.isRevToAlignment ) {
                        char *resSeqTmp = getNuclRevFragment(resSeq, res.dbLen, (NucleotideMatrix *) subMat);
                        ryId = getRYSeqId(res, querySeq, resSeqTmp, ryMap);
                        delete[] resSeqTmp;
                    }
                    else {
                        ryId = getRYSeqId(res, querySeq, resSeq, ryMap);
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
                    if (righty.isRevToAlignment) {
                        // Convert the reversed fragment to std::string
                        char *rightExtCandiSeqTmp = getNuclRevFragment(rightExtCandiSequ, rightExtCandiLen, (NucleotideMatrix *)subMat);
                        rightExtCandiSeq = std::string(rightExtCandiSeqTmp, rightExtCandiLen);
                        delete[] rightExtCandiSeqTmp;
                        // Here we calculate the epected similarity of the extension
                        isBadCandi = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNucRev, seqErrMatch, true);
                    } else {
                        // Directly convert the raw sequence to std::string
                        rightExtCandiSeq = std::string(rightExtCandiSequ, rightExtCandiLen);
                        // Here we calculate the epected similarity of the extension
                        isBadCandi = calcLikelihoodCorrection(rightLongestExt, righty, rightLongestExtSeq, rightExtCandiSeq, subDeamDiNuc, seqErrMatch, true);
                    }

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
                    if (lefty.isRevToAlignment) {
                        // Convert the reversed fragment to std::string
                        char *leftExtCandiSeqTmp = getNuclRevFragment(leftExtCandiSequ, leftExtCandiLen, (NucleotideMatrix *)subMat);
                        leftExtCandiSeq = std::string(leftExtCandiSeqTmp, leftExtCandiLen);
                        delete[] leftExtCandiSeqTmp;
                        // Here we calculate the epected similarity of the extension
                        isBadCandi = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNucRev, seqErrMatch, false);
                    } else {
                        // Directly convert the raw sequence to std::string
                        leftExtCandiSeq = std::string(leftExtCandiSequ, leftExtCandiLen);
                        // Here we calculate the epected similarity of the extension
                        isBadCandi = calcLikelihoodCorrection(leftLongestExt, lefty, leftLongestExtSeq, leftExtCandiSeq, subDeamDiNuc, seqErrMatch, false);
                    }

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

            mgeCandiRight.clear();
            mgeCandiLeft.clear();
// MGE marker end

            // if ( queryKey == 4631784 ){
            //     std::cerr << "mge found left?:\t" << mgeFoundLeft << std::endl;
            //     std::cerr << "mge found right?:\t" << mgeFoundRight << std::endl;
            // }


            // choose only reads
            std::vector<Matcher::result_t> reads;

            //DEBUG CORRECTION
            if ( queryKey == 4631784 ){
                std::cerr << ">" << "QUERY_4631784" << std::endl;
                std::cerr << querySeq << std::endl;
            }
            // DEBUG CORRECTION END

            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

                Matcher::result_t targetRead = alignments[alnIdx];     
                unsigned int tId = sequenceDbr->getId(targetRead.dbKey);
                char *tSeq = sequenceDbr->getData(tId, thread_idx);
                bool deleteTargetSeq = false;

                // bool isContig = sequenceDbr->getExtData(tId);

                // if ( isContig ){
                //     continue;
                // }

                if (targetRead.isRevToAlignment){
                    tSeq = getNuclRevFragment(tSeq, targetRead.dbLen, (NucleotideMatrix *) subMat);
                    deleteTargetSeq = true;
                }
                
                // float insideSeqId = getSubSeqId(targetRead, querySeq, tSeq);
                // float insideSeqThr = 0.95;

                float ryId = getRYSeqId(targetRead, querySeq,  tSeq, ryMap);
                targetRead.rySeqId = ryId;

                float rymerThresh = par.correctionThreshold;
                //float rymerThresh = 0.95;
                if ( targetRead.alnLength <= 100){
                    rymerThresh = (static_cast<float>(targetRead.alnLength) - 1) / static_cast<float>(targetRead.alnLength);
                    rymerThresh = std::floor(rymerThresh * 1000) / 1000;
                }


#ifdef DEBUGCORR
                if ( queryKey == 4631784 ){
                // DEBUG CORRECTION
                    #pragma omp critical
                    {
                        std::string targetOverlap;
                        for ( int i = targetRead.dbStartPos; i <= targetRead.dbEndPos; i++ )
                        {
                            targetOverlap += tSeq[i];
                        }
                        std::string targetAligned(querySeqLen, '-');
                        // Insert the overlapping region
                        for (int i = targetRead.qStartPos, k = 0; i <= targetRead.qEndPos; i++, k++) {
                            targetAligned[i] = targetOverlap[k];
                        }
                        // Print the target aligned sequence with an appropriate label
                        //std::cerr << targetAligned << "\t" << targetRead.dbKey << std::endl;

                        // if ( queryKey == 4631784 ){
                        //     std::cerr << ">" << targetRead.dbKey << "_" << tId << std::endl;
                        //     std::cerr << tSeq << std::endl;
                        // }

                        // Append the loop output to outputLine
                        targetAligned += "\t";
                        targetAligned += std::to_string(targetRead.dbKey);
                        targetAligned += "\t";
                        targetAligned += std::to_string(targetRead.seqId);
                        targetAligned += "\t";
                        // targetAligned += std::to_string(insideSeqId);
                        // targetAligned += "\t";
                        targetAligned += std::to_string(targetRead.rySeqId);
                        targetAligned += "\t";
                        targetAligned += std::to_string(rymerThresh);
                        targetAligned += "\t";
                        targetAligned += std::to_string(mgeFoundRight);
                        targetAligned += "\t";
                        targetAligned += std::to_string(mgeFoundLeft);
                        targetAligned += "\t";
                        targetAligned += std::to_string(alignments.size());
                        targetAligned += "\n";

                        if (outputFileOver.is_open()) {
                            outputFileOver << targetAligned;
                        }
                    }
                }
#endif

                if (deleteTargetSeq) {
                    delete[] tSeq;
                }

                if ( targetRead.rySeqId >= rymerThresh && !mgeFoundRight && targetRead.dbStartPos == 0 && static_cast<unsigned int>(targetRead.qEndPos) == (querySeqLen - 1)) {
                    // right extension
                    reads.push_back(alignments[alnIdx]);
                }
                else if ( targetRead.rySeqId >= rymerThresh && !mgeFoundLeft && alignments[alnIdx].qStartPos == 0 && static_cast<unsigned int>(alignments[alnIdx].dbEndPos) == (alignments[alnIdx].dbLen - 1)) {
                    // left extension
                    reads.push_back(alignments[alnIdx]);
                }
                else if (targetRead.rySeqId >= rymerThresh ){
                    // this is for the targets that align within the query and do not potentially extend anything, therefore the MGE finder does not make sense here
                    reads.push_back(alignments[alnIdx]);
                } 
            }
            alignments.clear();


            // Iterating through all alignments per query
            for ( size_t targ = 0; targ < reads.size(); targ++){   
                
                Matcher::result_t target = reads[targ];
                
                unsigned int targetId = sequenceDbr->getId(target.dbKey);
                char *tarSeq = sequenceDbr->getData(targetId, thread_idx);
                char *tSeq;

                //const bool targetWasExt = sequenceDbr->getExtData(targetId);

                // Get target sequence here to avoid passing the sequenceDbr
                unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId);
                bool deleteTargetSeq = false;

                if (useReverse[targetId]) {
                    tSeq = getNuclRevFragment(tarSeq, targetSeqLen, (NucleotideMatrix *) subMat);
                    deleteTargetSeq = true;
                }
                else {
                    tSeq = tarSeq;
                }

                float ryId = getRYSeqId(target, querySeq,  tSeq, ryMap);
                target.rySeqId = ryId;

                float rymerThresh = par.correctionThreshold;
                //float rymerThresh = 0.95;
                if ( target.alnLength <= 100){
                    rymerThresh = (static_cast<float>(target.alnLength) - 1) / static_cast<float>(target.alnLength);
                    rymerThresh = std::floor(rymerThresh * 1000) / 1000;
                }

                //if ( targetWasExt == false && isNotIdentity && target.rySeqId >= rymerThresh && target.seqId >= par.seqIdThr && target.alnLength >= 30 && ){
                //if ( targetWasExt == false && target.rySeqId >= rymerThresh && target.seqId >= par.seqIdThr && target.alnLength >= 30 ) {
                if ( target.rySeqId >= rymerThresh && target.seqId >= par.seqIdThr && target.alnLength >= 30 ) {

                    std::vector<int> subdeam_index(target.dbLen);
                    // Create lookup vector
                    for (size_t i = 0; i < 5; ++i) {
                        subdeam_index[i] = i;
                    }
                    for (size_t i = 5; i < target.dbLen - 5; ++i) {
                        subdeam_index[i] = 5;
                    }
                    for (size_t i = 0; i < 5; ++i) {
                        subdeam_index[target.dbLen - 5 + i] = 6 + i;
                    }

                    for ( unsigned int pos = 0; pos < target.alnLength; pos++ ){

                        int posInQuery = target.qStartPos + pos;

                        int targetBase = nucleotideMap[tSeq[target.dbStartPos + pos]];

                        queryCov[posInQuery][targetBase] += 1;

                        totalCov[posInQuery] += 1;
                        int deamT = subdeam_index[target.dbStartPos + pos];
                        deamVec[posInQuery].count[targetBase][deamT] += 1;
                        revCount[posInQuery].count[targetBase][deamT] += target.isRevToAlignment;
                    }
                }

                if ( deleteTargetSeq )
                    {
                        delete[] tSeq;
                    }
            }

            char * corrQuery = new char[querySeqLen];

            // For the Query Sequence 
            for ( unsigned int qPos = 0; qPos < querySeqLen; qPos++ )
            {
                int qBase = nucleotideMap[querySeq[qPos]];
                // Getting coverage of each base A,C,G,T at position "pos" in alignment (=overlap)
                if ( qWasExtended == false ){
                    queryCov[qPos][qBase] += 1;
                    totalCov[qPos] += 1;

                    if ( qPos < 5 ){
                        deamVec[qPos].count[qBase][qPos] += 1;
                    }
                    else if ( qPos >= querySeqLen - 5 ){
                        int queryIdx = querySeqLen - 1 - qPos;
                        deamVec[qPos].count[qBase][ 10 - queryIdx] += 1;
                    }
                    else{
                        deamVec[qPos].count[qBase][5] += 1;
                    }
                }
               
                if ( totalCov[qPos] <= 1){
                    corrQuery[qPos] = querySeq[qPos];
                }
                // if the sequence is a contig and has not been extended or corrected then just go by coverage:
                // else if ( qWasExtended == false ){
                //     unsigned int maxVal = 0;
                //     bool flags[4] = {false, false, false, false};
                //     size_t maxCount = 0;
                //     size_t maxIdx = qBase;

                //     // Find the max value in the current row
                //     for (size_t j = 0; j < queryCov[qPos].size(); ++j) {
                //         if (queryCov[qPos][j] > maxVal) {
                //             maxVal = queryCov[qPos][j];
                //         }
                //     }

                //     // Set flags for indices where the max value occurs
                //     for (size_t j = 0; j < queryCov[qPos].size(); ++j) {
                //         if (queryCov[qPos][j] == maxVal) {
                //             flags[j] = true;
                //             maxCount++;
                //         }
                //     }
                //     if ( (flags[qBase] && maxCount > 1) || (!flags[qBase] && maxCount > 1) ){
                //         corrQuery[qPos] = querySeq[qPos];
                //     }
                //     else{
                //         for (size_t j = 0; j < queryCov[qPos].size(); ++j) {
                //             if (queryCov[qPos][j] == maxVal) {
                //                 maxIdx = j;
                //                 break;
                //             }
                //         }
                //         corrQuery[qPos] = "ACGT"[maxIdx];
                //     }
                // }
                else {
                    int newBase = qBase;
                    int newBaseCandidate = mostLikeliBaseRead(qBase, qPos, deamVec, revCount, subDeamDiNuc, subDeamDiNucRev, seqErrMatch);
                    newBase = newBaseCandidate;

                    corrQuery[qPos] = "ACGT"[newBase];
                }
            }

            std::string corrStr(corrQuery, querySeqLen);
            corrStr.push_back('\n');
            delete[] corrQuery;


#ifdef DEBUG_CORR
        #pragma omp critical
        {
            std::string outputLine = std::to_string(queryKey) + "\t";

            std::ostringstream oss;
            for (unsigned int i = 0; i < totalCov.size(); i++) {
                oss << totalCov[i];
                if (i < totalCov.size() - 1) {
                    oss << " ";
                }
            }

            // Append the loop output to outputLine
            outputLine += oss.str();
            outputLine += "\n";

            if (outputFileAll.is_open()) {
                outputFileAll << outputLine;
            }
        }
#endif

            //if ( queryKey == 1250638 ){
            if ( queryKey == 4631784 ){
            //unsigned int max_cov = *std::max_element(totalCov.begin(), totalCov.end());
            //if ( max_cov >= 4 ){
            //if ( queryKey == 598571 ){
                std::cerr << ">Original\n" << querySeq;
                std::cerr << ">Corrected\n" << corrStr;
                std::cerr << "TotalCov:\n";
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << totalCov[i] << " ";
                }
                std::cerr << "\n";

                std::cerr << "QueryCov:\n";
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][0] << " ";
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][1] << " ";
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][2] << " ";
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][3] << " ";
                }
                std::cerr << "\n" << std::endl;
            }



           //resultWriter.writeData(corrQuery, querySeqLen, queryKey, thread_idx, qWasExtended);
            resultWriter.writeData(corrStr.c_str(), corrStr.size(), queryKey, thread_idx, qWasExtended);
            //corrQuery.clear();
            totalCov.clear();
            queryCov.clear();
            deamVec.clear();
            
        }
    } // end parallel

    // cleanup
    resultWriter.close(true);
    alnReader->close();
    delete [] wasExtended;
    delete alnReader;
    sequenceDbr->close();
    delete sequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    std::cout << "Correction end" << std::endl;
    return EXIT_SUCCESS;
}

int correction(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Do correction.\n";
    return doCorrection(par);
}

