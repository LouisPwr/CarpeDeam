#include "nuclassembleUtil.h"
const double SMOOTHING_VALUE = 0.001;

//#define DEBUGCORR
//#define DEBUCORR2

int mostLikeliBaseRead(const int baseInQuery, const unsigned int qIter, const std::vector<countDeamCov> & deamVec, const std::vector<countDeamCov> & countRevs, const std::vector<diNucleotideProb> & subDeamDiNuc, const std::vector<diNucleotideProb> & subDeamDiNucRev, const diNucleotideProb & seqErrMatch, bool wasCorr, unsigned int querySeqLen)
{
    // Initial likelihod of qBase (either A,C,G or T)
    //std::vector<float> baseFreqs = { 0.23554, 0.26446, 0.26446, 0.23554 };
    std::vector<float> baseFreqs = { 0.25, 0.25, 0.25, 0.25 };
    std::vector<long double> baseLikelis(4, 0.0); 

    // Precompute invariant logs
    std::vector<double> logBaseFreqs(4);
    std::vector<double> logQBaseErr(4);
    std::vector<double> logTBaseErr(4);

    unsigned int covA = 0;
    unsigned int covC = 0;
    unsigned int covG = 0;
    unsigned int covT = 0;
    for (size_t a = 0; a < subDeamDiNuc.size(); a++) {
        covA += deamVec[qIter].count[0][a];
    }
    for (size_t c = 0; c < subDeamDiNuc.size(); c++) {
        covC += deamVec[qIter].count[1][c];
    }
    for (size_t g = 0; g < subDeamDiNuc.size(); g++) {
        covG += deamVec[qIter].count[2][g];
    }
    for (size_t t = 0; t < subDeamDiNuc.size(); t++) {
        covT += deamVec[qIter].count[3][t];
    }
    std::vector<unsigned int> baseCovs = {covA, covC, covG, covT}; 


    double ctRatio = static_cast<double>(baseCovs[3]) /  (baseCovs[1] + baseCovs[3] + baseCovs[0] + baseCovs[2]); // Ensure floating-point division
    double gaRatio = static_cast<double>(baseCovs[0]) /  (baseCovs[1] + baseCovs[3] + baseCovs[0] + baseCovs[2]); // Ensure floating-point division

// maybe use
// double ctSum = baseCovs[1] + baseCovs[3];
// double gaSum = baseCovs[0] + baseCovs[2];

// double ctRatio = (ctSum != 0) ? static_cast<double>(baseCovs[3]) / ctSum  : 0.0;
// double gaRatio = (gaSum != 0) ? static_cast<double>(baseCovs[0]) / gaSum : 0.0;
 
    for (unsigned int q = 0; q < 4; q++) {
        logBaseFreqs[q] = std::log(baseFreqs[q]);  // Assuming baseFreqs are all 0.25
        //logQBaseErr[q] = std::log(seqErrMatch.p[q][baseInQuery]);
        logTBaseErr[q] = std::log(seqErrMatch.p[q][baseInQuery]);
        if ( wasCorr == true){             
            logQBaseErr[q] = std::log(seqErrMatch.p[q][baseInQuery]);
        }        
        else{
            if (ctRatio >= 0.4 || gaRatio >= 0.4 ) {
                //auto maxElementIt = std::max_element(baseCovs.begin(), baseCovs.end());
                //int maxIndex = std::distance(baseCovs.begin(), maxElementIt);
                return baseInQuery;
            } 
            if ( qIter < 5 ){
                double deamination = subDeamDiNuc[qIter].p[q][baseInQuery];
                double QBaseErr = std::max(deamination, SMOOTHING_VALUE);
                logQBaseErr[q] = std::log(QBaseErr);
            }
            else if ( qIter >= querySeqLen - 5 ){    
                double deamination = subDeamDiNuc[subDeamDiNuc.size() - (querySeqLen - qIter)].p[q][baseInQuery];
                double QBaseErr = std::max(deamination, SMOOTHING_VALUE);
                logQBaseErr[q] = std::log(QBaseErr);
            }
            else{
                double deamination = subDeamDiNuc[5].p[q][baseInQuery];
                double QBaseErr = std::max(deamination, SMOOTHING_VALUE);
                logQBaseErr[q] = std::log(QBaseErr);
            }
        }
    }

    // calculating the most likeli base
    for ( int qBase=0; qBase<4; qBase++ )
    {
        long double qBaseLik = 0;
        for (int tBase=0; tBase < 4; tBase++)
        {
            if ( baseCovs[tBase] == 0 ){
                //baseLikelis[tBase] = -std::numeric_limits<long double>::infinity(); // log(0)
                continue;
            }
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

                // int covDeam = deamVec[qIter].count[tBase][l];
                // int numReverse = countRevs[qIter].count[tBase][l];
                int covDeam = deamVec[qIter].count[tBase][l];
                int numReverse = countRevs[qIter].count[tBase][l];

                double logDeamPattern = std::log(deamPattern);
                double logDeamPatternRev = std::log(deamPatternRev);

                if (covDeam != 0) {
                    qBaseLik += (covDeam - numReverse) * (logTBaseErr[tBase] + logQBaseErr[qBase] + logDeamPattern);
                    qBaseLik += numReverse * (logTBaseErr[tBase] +  logQBaseErr[qBase] + logDeamPatternRev);
                }
                //qBaseLik += (covDeam - numReverse) * (logTBaseErr[tBase] + logQBaseErr[qBase] + logDeamPattern); 
            }
        }
        baseLikelis[qBase] = qBaseLik;
    }

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
            //auto startQuery = std::chrono::high_resolution_clock::now();
            // double sumLikeli = 0;

            progress.updateProgress();

            unsigned int queryKey = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id, thread_idx);
            unsigned int querySeqLen = sequenceDbr->getSeqLen(id);
            //std::string query(querySeq, querySeqLen); // no /n/0

            char *alnData = alnReader->getDataByDBKey(queryKey, thread_idx);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);

            bool qWasExtended = sequenceDbr->getExtData(id);

            float avCov = 0;

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
                avCov += alignments[alnIdx].alnLength ;
            }
            avCov = static_cast<float>(avCov)/querySeqLen;

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

            // choose only reads
            std::vector<Matcher::result_t> reads;
            reads.reserve(300);

            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

                Matcher::result_t targetRead = alignments[alnIdx];     
                unsigned int tId = sequenceDbr->getId(targetRead.dbKey);
                char *tSeq = sequenceDbr->getData(tId, thread_idx);
                bool deleteTargetSeq = false;

                bool isContig = sequenceDbr->getExtData(tId);

                if ( isContig ){
                    continue;
                }

                if (targetRead.isRevToAlignment){
                    tSeq = getNuclRevFragment(tSeq, targetRead.dbLen, (NucleotideMatrix *) subMat);
                    deleteTargetSeq = true;
                }
                
                // float insideSeqId = getSubSeqId(targetRead, querySeq, tSeq);
                // float insideSeqThr = 0.95;

                float ryId = getRYSeqId(targetRead, querySeq,  tSeq, ryMap);
                targetRead.rySeqId = ryId;

                float rymerThresh = par.corrReadsRySeqId;
                if ( targetRead.alnLength <= 100){
                    rymerThresh = (static_cast<float>(targetRead.alnLength) - 1) / static_cast<float>(targetRead.alnLength);
                    rymerThresh = std::floor(rymerThresh * 1000) / 1000;
                }

                if (deleteTargetSeq) {
                    delete[] tSeq;
                }

                //std::pair<bool, bool> mgefound = mgeFinder(alignments, sequenceDbr, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat, subDeamDiNuc, subDeamDiNucRev, seqErrMatch);
                std::pair<bool,bool> mgefound(false, false);


                if ( targetRead.rySeqId >= rymerThresh && !mgefound.second && targetRead.dbStartPos == 0 && static_cast<unsigned int>(targetRead.qEndPos) == (querySeqLen - 1)) {
                    // right extension
                    reads.push_back(alignments[alnIdx]);
                }
                else if ( targetRead.rySeqId >= rymerThresh && !mgefound.first && targetRead.qStartPos == 0 && static_cast<unsigned int>(targetRead.dbEndPos) == (targetRead.dbLen - 1)) {
                    // left extension
                    reads.push_back(alignments[alnIdx]);
                }
                else if (targetRead.rySeqId >= rymerThresh && avCov < 50 ){
                    // this is for the targets that align within the query and do not potentially extend anything, therefore the MGE finder does not make sense here
                    reads.push_back(alignments[alnIdx]);
                } 
            }
            //alignments.clear();


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

                float rymerThresh = par.corrReadsRySeqId;
                if ( target.alnLength <= 100){
                    rymerThresh = (static_cast<float>(target.alnLength) - 1) / static_cast<float>(target.alnLength);
                    rymerThresh = std::floor(rymerThresh * 1000) / 1000;
                }

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

                // if ( qWasExtended == false ){
                    // queryCov[qPos][qBase] += 1;
                    // totalCov[qPos] += 1;

                    // if ( qPos < 5 ){
                    //     deamVec[qPos].count[qBase][qPos] += 1;
                    // }
                    // else if ( qPos >= querySeqLen - 5 ){
                    //     int queryIdx = querySeqLen - 1 - qPos;
                    //     deamVec[qPos].count[qBase][ 10 - queryIdx] += 1;
                    // }
                    // else{
                    //     deamVec[qPos].count[qBase][5] += 1;
                    // }
                // }
               
                if ( totalCov[qPos] <= 1 ){
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
                    //auto startLikeli = std::chrono::high_resolution_clock::now();
                    int newBase = qBase;
                    int newBaseCandidate = mostLikeliBaseRead(qBase, qPos, deamVec, revCount, subDeamDiNuc, subDeamDiNucRev, seqErrMatch, qWasExtended, querySeqLen);
                    newBase = newBaseCandidate;

                    corrQuery[qPos] = "ACGT"[newBase];
                }
            }

            std::string corrStr(corrQuery, querySeqLen);
            corrStr.push_back('\n');
            delete[] corrQuery;

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

