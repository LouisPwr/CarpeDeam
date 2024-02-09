#include "nuclassembleUtil.h"
const double SMOOTHING_VALUE = 0.0000001;


int mostLikeliBaseRead(int qIter, diNucleotideProb & match, diNucleotideProb & mismatch, std::vector<countDeamCov> & deamVec, std::vector<countDeamCov> & countRevs, std::vector<diNucleotideProb> & subDeamDiNuc, std::vector<diNucleotideProb> & subDeamDiNucRev)
{
    std::vector<double> baseLikelis; 

    // Initial likelihod of qBase (either A,C,G or T)
    std::vector<float> baseFreqs = { 0.23554, 0.26446, 0.26446, 0.23554 };

    for ( int qBase=0; qBase<4; qBase++ )
    {
        double qBaseLik = 0;
        //double qBaseLogLik = 0;
        for (int tBase=0; tBase < 4; tBase++)
        {
            double seqLik = (qBase == tBase) ? match.p[qBase][tBase] : mismatch.p[qBase][tBase];

            // The vector deamVec contains the counts of deaminations per position;
            // Iterating through it to make the amount of multiplications we need to do. 
            for ( unsigned int l = 0; l < subDeamDiNuc.size(); l++)
            {
                double deamPattern = subDeamDiNuc[l].p[qBase][tBase];
                double deamPatternRev = subDeamDiNucRev[l].p[qBase][tBase];
                // this can only be true in case of a mismatch of Pur<->Pyr
                deamPattern = (deamPattern == 0.0) ? SMOOTHING_VALUE : deamPattern;
                deamPatternRev = (deamPatternRev == 0.0) ? SMOOTHING_VALUE : deamPatternRev;

                int covDeam = deamVec[qIter].count[tBase][l];
                int numReverse = countRevs[qIter].count[tBase][l];
                qBaseLik += (covDeam - numReverse)*( log( baseFreqs[qBase] ) + log(deamPattern) + log(seqLik) );
                qBaseLik += numReverse*( log( baseFreqs[qBase] ) + log(deamPatternRev) + log(seqLik) );
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

/* int mostLikeliBaseContig(int contigBase, int qIter, diNucleotideProb & match, diNucleotideProb & mismatch, std::vector<countDeamCov> & deamVec, std::vector<countDeamCov> & countRevs, std::vector<diNucleotideProb> & subDeamDiNuc, std::vector<diNucleotideProb> & subDeamDiNucRev)
{
    std::vector<double> baseLikelis; 

    // Initial likelihod of qBase (either A,C,G or T)
    std::vector<float> baseFreqs = { 0.23554, 0.26446, 0.26446, 0.23554 };
    for ( int qBase=0; qBase<4; qBase++ ){
        double qBaseLik = 0;
        double seqLik = match.p[qBase][contigBase];
        for (int tBase=0; tBase < 4; tBase++){
            // The vector deamVec contains the counts of deaminations per position;
            // Iterating through it to make the amount of multiplications we need to do. 
            for ( unsigned int l = 0; l < subDeamDiNuc.size(); l++){
                double deamPattern = subDeamDiNuc[l].p[qBase][tBase];
                double deamPatternRev = subDeamDiNucRev[l].p[qBase][tBase];
                // this can only be true in case of a mismatch of Pur<->Pyr
                if ( deamPattern == 0.0 )
                {
                    deamPattern = 1.0;
                }
                if ( deamPatternRev == 0.0 )
                {
                    deamPatternRev = 1.0;
                }

                int covDeam = deamVec[qIter].count[tBase][l];
                int numReverse = countRevs[qIter].count[tBase][l];
                qBaseLik += (covDeam - numReverse)*( log( baseFreqs[qBase] ) + log(deamPattern) + log(seqLik) );
                qBaseLik += numReverse*( log( baseFreqs[qBase] ) + log(deamPatternRev) + log(seqLik) );
            }
        }
        baseLikelis.push_back(qBaseLik);
    }

    auto maxLik = std::max_element(baseLikelis.begin(), baseLikelis.end());
    int idxMax = std::distance(baseLikelis.begin(), maxLik);

    return idxMax;
} */



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
    // if ( userInput.empty() ){
    //     createNoDamageMatrix("none.prof");
    //     high5 = "none.prof";
    //     high3 = "none.prof";
    // }

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


        // define length here
        //int maxFragLen = 250; 
        //std::vector<std::vector<diNucleotideProb>> allDeam(maxFragLen);
        std::vector<diNucleotideProb> subDeamDiNuc(11);
        //Substitution rates due to deamination
        std::vector<substitutionRates> sub5p;
        std::vector<substitutionRates> sub3p;

        std::vector<diNucleotideProb> subDeamDiNucRev;

        //std::cerr << "Before initDeamProbabilities" << std::endl;
        initDeamProbabilities(high5, high3, sub5p, sub3p, subDeamDiNuc, subDeamDiNucRev);

        //std::vector<diNucleotideProb> subDeamDiNucRev;

        diNucleotideProb seqErrMatch;
        diNucleotideProb seqErrMis;

        long double seqErrCorrection= 0.01;
        getSeqErrorProf(seqErrMatch, seqErrMis, seqErrCorrection);

        float rymerThresh = par.correctionThreshold;


        int reverse = 0;
        int good = 0;

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

            bool qWasExtended = false;

            ////std::cerr << "query id: " << queryKey << std::endl;
            //std::cerr << "\n" << "new query: " << id << std::endl;
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
                        reverse += 1;

                    } else {
                        useReverse[sequenceDbr->getId(alignments[alnIdx].dbKey)] = false;
                        alignments[alnIdx].isRevToAlignment = false;
                        good += 1;
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

            qWasExtended = sequenceDbr->getExtData(id);

            // Iterating through all alignments per query
            for ( size_t targ = 0; targ < alignments.size(); targ++){   
                
                Matcher::result_t target = alignments[targ];

                //const bool notInside = target.dbLen != target.alnLength;
                //const bool rightStart = target.dbStartPos == 0;
                //const bool leftStart = target.qStartPos == 0;

                unsigned int targetId = sequenceDbr->getId(target.dbKey);
                char *tarSeq = sequenceDbr->getData(targetId, thread_idx);
                char *tSeq;

                //const bool notInside = true;
                //const bool rightStart = true;
                //const bool leftStart = true;
                const bool isNotIdentity = (target.dbKey != queryKey);

                const bool targetWasExt = sequenceDbr->getExtData(targetId);

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
        //}

                float ryId = getRYSeqId(target, querySeq,  tSeq, ryMap);
                target.rySeqId = ryId;

                //float subSeqId = getSubSeqId(target, querySeq, tSeq);

                //std::cerr << "RySeqId:\t" << target.rySeqId;
                //std::cerr << "\n" << std::endl;

                if ( targetWasExt == false && isNotIdentity && target.rySeqId >= rymerThresh && target.alnLength >= 15 ){

                    //std::cerr << ">" << target.dbKey << std::endl;
                    //std::cerr << tSeq; 

                    for ( int pos = target.qStartPos; pos <= target.qEndPos; pos++ ){

                        int queryIter = pos - target.qStartPos;
                        int targetBase = nucleotideMap[tSeq[target.dbStartPos + queryIter]];
                        // Getting coverage of each base A,C,G,T at position "pos" in alignment (=overlap)
                        queryCov[pos][targetBase] += 1;
                        totalCov[pos] += 1;

                        int relativePos = target.dbStartPos + queryIter;
                        if ( relativePos < 5 ){
                            deamVec[pos].count[targetBase][relativePos] += 1;
                            revCount[pos].count[targetBase][relativePos] += target.isRevToAlignment;
                            //tProbs = subDeamDiNuc[relativePos];
                        }
                        else if ( relativePos > target.dbEndPos - 5 ){
                            deamVec[pos].count[targetBase][ 10 - ( target.dbEndPos - relativePos )] += 1;
                            revCount[pos].count[targetBase][ 10 - ( target.dbEndPos - relativePos )] += target.isRevToAlignment;
                            //tProbs = subDeamDiNuc[ 10 - (target.dbEndPos - relativePos) ];
                        }
                        else{
                            deamVec[pos].count[targetBase][5] += 1;
                            revCount[pos].count[targetBase][5] += target.isRevToAlignment;
                            //tProbs = subDeamDiNuc[5];
                        }
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
                // If query is a read, add it to the coverage as it should be included in the voting
                if ( totalCov[qPos] > 0 && qWasExtended == false )
                {
                    // Getting coverage of each base A,C,G,T at position "pos" in alignment (=overlap)
                    queryCov[qPos][qBase] += 1;
                    totalCov[qPos] += 1;

                    if ( qPos < 5 ){
                        deamVec[qPos].count[qBase][qPos] += 1;
                    }
                    else if ( qPos >= querySeqLen - 5 ){
                        deamVec[qPos].count[qBase][ 10 - ( (querySeqLen-1) - qPos )] += 1;
                    }
                    else{
                        deamVec[qPos].count[qBase][5] += 1;
                    }
                }

                if ( totalCov[qPos] <= 1){
                    corrQuery[qPos] = querySeq[qPos];
                }
                else {
                    int newBase = mostLikeliBaseRead(qPos, seqErrMatch, seqErrMis, deamVec, revCount, subDeamDiNuc, subDeamDiNucRev);
                    corrQuery[qPos] = "ACGT"[newBase];
                    // if ( newBase != qBase )
                    // {
                    //     std::cerr << "Original:\t" << querySeq << "qPos:\t" << qPos << "\tquerySeqLen:\t" << querySeqLen << std::endl;
                    //     std::cerr << "A:" << queryCov[qPos][0] << " C:" << queryCov[qPos][1] << " G:" << queryCov[qPos][2] << " T:" << queryCov[qPos][3] << std::endl;
                    //     std::cerr << "Total Cov: " << totalCov[qPos] << std::endl;
                    //     std::cerr << "old:\t" << qBase << "\tnew:\t" << newBase << std::endl;
                    // }
                }
            }

            std::string corrStr(corrQuery, querySeqLen);
            corrStr.push_back('\n');
            delete[] corrQuery;

            if ( false )
            {
                std::cerr << ">Original\n" << querySeq;
                std::cerr << ">Corrected\n" << corrStr;
                std::cerr << "TotalCov:\n";
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    //std::cerr << totalCov[i];
                }
                std::cerr << "\n";

                std::cerr << "QueryCov:\n";
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][0];
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][1];
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][2];
                }
                std::cerr << std::endl;
                for (unsigned int i = 0; i < totalCov.size(); i++)
                {
                    std::cerr << queryCov[i][3];
                }
                std::cerr << "\n" << std::endl;
            }



           //resultWriter.writeData(corrQuery, querySeqLen, queryKey, thread_idx, qWasExtended);
            resultWriter.writeData(corrStr.c_str(), corrStr.size(), queryKey, thread_idx, qWasExtended);
            //corrQuery.clear();
            totalCov.clear();
            queryCov.clear();
            deamVec.clear();
            

/*             std::cerr << "rev: " << reverse << std::endl;
            std::cerr << "good: " << good << std::endl; */

            // Need to change the deams as they change for each target !!! //
            /////////////////////////////////////////////////////////////////////////////////////////////////////////
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

/* int helloworld()
{
    std::cout << "hello world" << std::endl;
    exit(1);
    return 0;
} */

int correction(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Do correction.\n";
    return doCorrection(par);
}

