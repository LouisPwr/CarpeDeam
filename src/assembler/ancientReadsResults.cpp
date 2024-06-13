#include "nuclassembleUtil.h"

//#define DEBUGEXT

//#define LIKELI

std::vector<unsigned int> getMaxAlnLen(std::vector<Matcher::result_t> &alignments, unsigned int & queryKey)
{
    unsigned int rightLen = 0;
    unsigned int leftLen = 0;
    std::vector<unsigned int> maxAlnLen;
    int left = 0;

    for (size_t r = 0; r < alignments.size(); r++)
    {

        //const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 &&  res.qStartPos == 0 );
        bool notInside = alignments[r].dbLen != alignments[r].alnLength;
        const bool rightStart = alignments[r].dbStartPos == 0;
        const bool leftStart = alignments[r].qStartPos == 0;
        const bool isNotIdentity = (alignments[r].dbKey != queryKey);

        if (rightStart && notInside && isNotIdentity && unsigned(alignments[r].qEndPos) == (alignments[r].qLen - 1)) {
            //right extension
            if (alignments[r].alnLength > rightLen)
            {
                rightLen = alignments[r].alnLength;
            }

        }
        else if (leftStart && notInside && isNotIdentity && unsigned(alignments[r].dbEndPos) == (alignments[r].dbLen - 1)) {
            //left extension
            left += 1;
            if (alignments[r].alnLength > leftLen)
            {
                leftLen = alignments[r].alnLength;
            }
        }
    }
    //std::cerr << "count left " << left << std::endl;

    maxAlnLen.push_back(leftLen);
    maxAlnLen.push_back(rightLen);

    return maxAlnLen;
}


scorePerRes r_s_pair(Matcher::result_t res, std::string & consensus, char* targetSeq, unsigned int querySeqLen, std::vector<diNucleotideProb> &subDeamDiNuc, std::vector<diNucleotideProb> &subDeamDiNucRev, unsigned int & maxLeft, unsigned int & maxRight, float randAlnPenal, diNucleotideProb & seqErrMatch, float excessPenal)
{
    scorePerRes fiveToThree;
    //scorePerRes threeToFive;

    fiveToThree.r = res;
    //threeToFive.r = res;

    unsigned int maxOverlap = maxRight;

    if (unsigned(res.qStartPos) == 0 && unsigned(res.dbEndPos) == (res.dbLen - 1)) {
        maxOverlap = maxLeft;
    }

    std::vector<diNucleotideProb> subDeamDiNucRef = res.isRevToAlignment ? subDeamDiNucRev : subDeamDiNuc;

    //calcLikelihood(fiveToThree, querySeq, targetSeq, subDeamDiNucRef, maxOverlap, randAlnPenal, seqErrMatch, excessPenal);
    calcLikelihoodConsensus(fiveToThree, consensus, querySeqLen, targetSeq, subDeamDiNucRef, maxOverlap, randAlnPenal, seqErrMatch, excessPenal);
    // Modifying ThreeToFive
    //std::cerr << " REV " << std::endl;
    //calc_likelihood(threeToFive, querySeq, targetSeq, subDeamDiNucRev, maxOverlap, isRightOverlap, randAlnPenal, seqErrMatch, seqErrMis);
    return fiveToThree;
}


// For extension with reads
typedef std::priority_queue<scorePerRes, std::vector<scorePerRes>, CompareNuclResultByScoreReads> QueueByScoreNuclReads;

Matcher::result_t selectNuclFragmentToExtendReads(QueueByScoreNuclReads &alignments, unsigned int queryKey) {
    // results are ordered by score
    while (alignments.empty() == false){
        Matcher::result_t res = alignments.top().r;
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 &&  res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);

        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}



int doNuclAssembly1(LocalParameters &par) {
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

    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);
    EvalueComputation evaluer(sequenceDbr->getAminoAcidDBSize(), subMat);

    unsigned char * wasExtended = new unsigned char[sequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+sequenceDbr->getSize(), 0);
    Debug::Progress progress(sequenceDbr->getSize());

    std::string userInput = par.ancientDamagePath; // This should come from the user

    std::string high5 = userInput + "5p.prof";
    std::string high3 = userInput + "3p.prof";

    //unsigned int totalRight = 0; 
    //unsigned int rightExtLongest = 0;

    //unsigned int totalLeft = 0; 
    //unsigned int leftExtLongest = 0;

    //unsigned int wasReverse = 0;
    //unsigned int totalQuer = 0;

#ifdef DEBUGEXT
    auto currentTime = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(currentTime);
    std::stringstream one;
    one << "likelihood_all_" << time << ".tsv";
    std::string filenameOne = one.str();

    std::stringstream two;
    two << "ext_reads_" << time << ".tsv";
    std::string filenameTwo = two.str();

    std::ofstream outputFileOne;  // Declare the ofstream object globally if possible
    outputFileOne.open(filenameOne, std::ios::app);  // Open once, append mode

    std::ofstream outputFileTwo;  // Declare the ofstream object globally if possible
    outputFileTwo.open(filenameTwo, std::ios::app);  // Open once, append mode

#endif

//#pragma omp parallel reduction(+:totalRight,totalLeft,rightExtLongest,leftExtLongest)
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

        long double seqErrCorrection= 0.001;
        getSeqErrorProf(seqErrMatch, seqErrCorrection);

        float randAlnPenal = par.randomAlignPenal;


#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            
            progress.updateProgress();

            unsigned int queryKey = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id, thread_idx);
            unsigned int querySeqLen = sequenceDbr->getSeqLen(id);
            std::string query(querySeq, querySeqLen); // no /n/0

            char *alnData = alnReader->getDataByDBKey(queryKey, thread_idx);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);

            bool queryCouldBeExtended = false;

            //CompareNuclResultByScore comparer(sequenceDbr, querySeq, deamMatrix, thread_idx, (NucleotideMatrix*)subMat);
            QueueByScoreNuclReads alnQueueReads;

            std::vector<scorePerRes> likScores;

            std::vector<Matcher::result_t> notContig;
            //std::vector<Matcher::result_t> contig;

            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {

                int rawScore = static_cast<int>(evaluer.computeRawScoreFromBitScore(alignments[alnIdx].score) + 0.5);
                float scorePerCol = static_cast<float>(rawScore) / static_cast<float>(alignments[alnIdx].alnLength + 0.5);

                //float alnLen = static_cast<float>(alignments[alnIdx].alnLength);
                //float ids = static_cast<float>(alignments[alnIdx].seqId) * alnLen;
                //alignments[alnIdx].seqId = ids / (alnLen + 0.5);
                alignments[alnIdx].score = static_cast<int>(scorePerCol*100);

                if (seqType == Parameters::DBTYPE_NUCLEOTIDES) {
                    if (alignments[alnIdx].qStartPos > alignments[alnIdx].qEndPos) {
                        useReverse[sequenceDbr->getId(alignments[alnIdx].dbKey)] = true;
                        std::swap(alignments[alnIdx].qStartPos, alignments[alnIdx].qEndPos);
                        unsigned int dbStartPos = alignments[alnIdx].dbStartPos;
                        alignments[alnIdx].dbStartPos = alignments[alnIdx].dbLen - alignments[alnIdx].dbEndPos - 1;
                        alignments[alnIdx].dbEndPos= alignments[alnIdx].dbLen - dbStartPos - 1;
                        alignments[alnIdx].isRevToAlignment = true;
                        // #pragma omp atomic
                        // wasReverse += 1;
                        // std::cerr << "Query:\t" << querySeq << std::endl;
                        // char *targetRev = sequenceDbr->getData(sequenceDbr->getId(alignments[alnIdx].dbKey), thread_idx);
                        // std::cerr << "Target:\t" << targetRev << std::endl;

                    } else {
                        useReverse[sequenceDbr->getId(alignments[alnIdx].dbKey)] = false;
                        alignments[alnIdx].isRevToAlignment = false;
                    }
                    // #pragma omp atomic
                    // totalQuer += 1;
                }

                if (alignments.size() > 1){
                    __sync_or_and_fetch(&wasExtended[sequenceDbr->getId(alignments[alnIdx].dbKey)], static_cast<unsigned char>(0x40));
                }
                                            
            }

            // update sequence identity
            for (unsigned int idx = 0; idx < alignments.size(); idx++){
                // now retrieve the other candidate extensions and get their sequences

                unsigned int aln2updateId = sequenceDbr->getId(alignments[idx].dbKey);
                unsigned int aln2updateLen = sequenceDbr->getSeqLen(aln2updateId);

                if ( aln2updateId == queryKey ){
                    continue;
                }

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

                // calculate the sequence identity
                int idCnt = 0;
                int idRyCnt = 0;
                for (int i = alignments[idx].qStartPos; i <= alignments[idx].qEndPos; i++) {
                    idCnt += (querySeq[i] == aln2updateSeq[alignments[idx].dbStartPos + (i-alignments[idx].qStartPos)]) ? 1 : 0;
                    idRyCnt += (ryMap[querySeq[i]] == ryMap[aln2updateSeq[alignments[idx].dbStartPos + (i-alignments[idx].qStartPos)]]) ? 1 : 0;
                }
                float seqId = static_cast<float>(idCnt) / alignments[idx].alnLength;
                float rySeqId = static_cast<float>(idRyCnt) / alignments[idx].alnLength;

                    // std::cerr << "alnLen " << alignments[idx].alnLength << std::endl;
                    // std::cerr << "idRyCnt " << idRyCnt << std::endl;
                    // std::cerr << "rySeqId " << rySeqId << std::endl;
                    // std::cerr << "idCnt " << idCnt << std::endl;
                    // std::cerr << "seqId " << seqId << std::endl;
                    // std::cerr << "query " << querySeq << std::endl;
                    // std::cerr << "aln2updateSequ " << aln2updateSequ << std::endl;
                    // std::cerr << "aln2updateSeq " << aln2updateSeq << std::endl;

                alignments[idx].seqId = seqId;
                alignments[idx].rySeqId = rySeqId;
            }

            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {
                bool noOffset = (alignments[alnIdx].dbLen - alignments[alnIdx].alnLength) == 0;
                unsigned int targetId = sequenceDbr->getId(alignments[alnIdx].dbKey);
                bool isContig = sequenceDbr->getExtData(targetId);
                //if ( isContig == false && alignments[alnIdx].alnLength >= 30 && alignments[alnIdx].seqId >= par.seqIdThr ){
                if ( isContig == false && alignments[alnIdx].alnLength >= 30 && alignments[alnIdx].seqId >= par.seqIdThr && !noOffset ){
                    notContig.push_back(alignments[alnIdx]);
                }
            }
            alignments.clear();

            std::string consensus = std::string(querySeq, querySeqLen);
            unsigned int maxAlnLeft = 0;
            unsigned int maxAlnRight = 0;
    
            consensus = consensusCaller(notContig, sequenceDbr, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat);
            updateSeqIdConsensusReads(notContig, sequenceDbr, consensus, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat, maxAlnLeft, maxAlnRight); 

            // ONLY DAMAGE AWARE EXTENSION

            // std::vector<unsigned int> maxAlign = getMaxAlnLen(notContig, queryKey);
            // unsigned int maxAlnLeft = maxAlign[0];
            // unsigned int maxAlnRight = maxAlign[1];
            // maxAlign.clear();

            for (size_t alnIdx = 0; alnIdx < notContig.size(); alnIdx++) {

                // Get target sequence here to avoid passing the sequenceDbr
                unsigned int targetId = sequenceDbr->getId(notContig[alnIdx].dbKey);
                char *targetSeq = sequenceDbr->getData(targetId, thread_idx);
                unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId);

                //calculate rySeqId for each result_t
                bool deleteTargetSeq = false; 
                if (notContig[alnIdx].isRevToAlignment)
                {
                    targetSeq = getNuclRevFragment(targetSeq, targetSeqLen, (NucleotideMatrix *) subMat);
                    deleteTargetSeq = true; 
                }

                // float ryId = getRYSeqId(notContig[alnIdx], querySeq,  targetSeq, ryMap);
                // ryId = std::round(ryId * 1000.0f) / 1000.0f;
                // notContig[alnIdx].rySeqId = ryId;

                bool notInside = notContig[alnIdx].dbLen != notContig[alnIdx].alnLength;
                const bool rightStart = notContig[alnIdx].dbStartPos == 0;
                const bool leftStart = notContig[alnIdx].qStartPos == 0;
                const bool isNotIdentity = (notContig[alnIdx].dbKey != queryKey);
                //const bool notRightStartAndLeftStart = !(notContig[alnIdx].dbStartPos == 0 && notContig[alnIdx].qStartPos == 0 );


#ifdef DEBUGEXT
                    std::string outputLine = std::to_string(queryKey) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(par.likelihoodThreshold) + "\t" +
                                                std::to_string(notContig[alnIdx].seqId) + "\t" +
                                                std::to_string(notContig[alnIdx].rySeqId) + "\t" +
                                                std::to_string(notInside) + "\t" +
                                                std::to_string(rightStart) + "\t" +
                                                std::to_string(leftStart) + "\t" +
                                                std::to_string(isNotIdentity) + "\t" +
                                                std::to_string(notRightStartAndLeftStart) + "\n";

                        #pragma omp critical
                        {
                            if (outputFileOne.is_open()) {
                                outputFileOne << outputLine;
                            }
                        }
#endif

                if ( (rightStart || leftStart) && notInside && isNotIdentity && notContig[alnIdx].rySeqId >= par.rySeqIdThr && notContig[alnIdx].seqId >= par.seqIdThr)
                {
                    // #pragma omp critical
                    // std::cerr << "Before FIRST" << std::endl;
                    // #pragma omp critical
                    // std::cerr << querySeq;
                    scorePerRes toAdd = r_s_pair(notContig[alnIdx], consensus, targetSeq, querySeqLen, subDeamDiNuc, subDeamDiNucRev, maxAlnLeft, maxAlnRight, randAlnPenal, seqErrMatch, par.excessPenal);
                    // #pragma omp critical
                    // std::cerr << "After FIRST" << std::endl;


#ifdef DEBUGEXT
                    std::string outputLine = std::to_string(toAdd.sRatio) + "\t" +
                                                std::to_string(toAdd.sLenNorm) + "\n";
                    #pragma omp critical
                    {
                        if (outputFileOne.is_open()) {
                            outputFileOne << outputLine;
                        }
                    }
#endif

                    if ( toAdd.sRatio > par.likelihoodThreshold ) {
                        alnQueueReads.push( toAdd );
                    }
                }

                if (deleteTargetSeq) {
                    delete[] targetSeq;
                }
            }

            std::vector<Matcher::result_t> tmpAlignmentsReads;
            tmpAlignmentsReads.reserve(notContig.size());
            // Louis was here
            //unsigned int alignment_counter = 0;

            while (!alnQueueReads.empty()) {

                unsigned int leftQueryOffset = 0;
                unsigned int rightQueryOffset = 0;
                tmpAlignmentsReads.clear();
                Matcher::result_t besttHitToExtend;

                while ((besttHitToExtend = selectNuclFragmentToExtendReads(alnQueueReads, queryKey)).dbKey != UINT_MAX) {

                    unsigned int targetId = sequenceDbr->getId(besttHitToExtend.dbKey);
                    if (targetId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Could not find targetId  " << besttHitToExtend.dbKey
                                            << " in database " << sequenceDbr->getDataFileName() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *targetSeq = sequenceDbr->getData(targetId, thread_idx);
                    unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId) ;

                    // check if alignment still make sense (can extend the query)
                    if (besttHitToExtend.dbStartPos == 0) {
                        if ((targetSeqLen - (besttHitToExtend.dbEndPos + 1)) <= rightQueryOffset) {
                            continue;
                        }
                    } else if (besttHitToExtend.qStartPos == 0) {
                        if (besttHitToExtend.dbStartPos <= static_cast<int>(leftQueryOffset)) {
                            continue;
                        }
                    }
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x10));

                    unsigned int dbStartPos = besttHitToExtend.dbStartPos;
                    unsigned int dbEndPos = besttHitToExtend.dbEndPos;
                    unsigned int qStartPos = besttHitToExtend.qStartPos;
                    unsigned int qEndPos = besttHitToExtend.qEndPos;
                    if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                        //right extension

                        if(rightQueryOffset > 0) {
                            tmpAlignmentsReads.push_back(besttHitToExtend);
                            continue;
                        }

                        unsigned int fragLen = targetSeqLen - (dbEndPos + 1);

                        if (query.size() + fragLen >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Ignore extension because of length limitation for sequence: " \
                                                  << queryKey << ". Max length allowed would be " << par.maxSeqLen << "\n";
                            break;
                        }

                        std::string fragment;
                        if (useReverse[targetId]) {
                            char *cfragment = getNuclRevFragment(targetSeq, fragLen, (NucleotideMatrix *) subMat);
                            fragment = std::string(cfragment, fragLen);
                            delete[] cfragment;
                        }
                        else
                           fragment = std::string(targetSeq + dbEndPos + 1, fragLen);

                        // #pragma omp atomic
                        // totalRight += 1;
                        // if ( maxAlnRight == besttHitToExtend.alnLength )
                        // {
                        //     #pragma omp atomic
                        //     rightExtLongest += 1;
                        // }
                        // std::cerr << query << std::endl;
                        // std::cerr << besttHitToExtend.seqId << std::endl;

                        query += fragment;
                        rightQueryOffset += fragLen;

#ifdef DEBUGEXT
                    std::string outputLine = std::to_string(queryKey) + "\t" +
                                                std::to_string(besttHitToExtend.alnLength) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(besttHitToExtend.dbLen) + "\t" +
                                                std::to_string(besttHitToExtend.seqId) + "\t" +
                                                std::to_string(besttHitToExtend.rySeqId) + "\t" +
                                                "1\t" +  // Assuming this is a constant value
                                                std::to_string(query.length()) + "\t\n";

                        #pragma omp critical
                        {
                            if (outputFileTwo.is_open()) {
                                outputFileTwo << outputLine;
                            }
                        }
#endif

                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));

                    }
                    else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                        //left extension

                        if(leftQueryOffset > 0) {
                            tmpAlignmentsReads.push_back(besttHitToExtend);
                            continue;
                        }

                        unsigned int fragLen = dbStartPos;
                        if (query.size() + fragLen >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Ignore extension because of length limitation for sequence: " \
                                                  << queryKey << ". Max length allowed would be " << par.maxSeqLen << "\n";
                            break;
                        }

                        std::string fragment;
                        if (useReverse[targetId]) {
                            char *cfragment = getNuclRevFragment(targetSeq + (targetSeqLen - dbStartPos), fragLen, (NucleotideMatrix *) subMat);
                            fragment = std::string(cfragment, fragLen);
                            delete[] cfragment;
                        }
                        else
                            fragment = std::string(targetSeq, fragLen);
                        
                        // #pragma omp atomic
                        // totalLeft += 1;
                        // if ( maxAlnLeft == besttHitToExtend.alnLength )
                        // {
                        //     #pragma omp atomic
                        //     leftExtLongest += 1;
                        // }

                        query = fragment + query;
                        leftQueryOffset += fragLen;

#ifdef DEBUGEXT
                    std::string outputLine = std::to_string(queryKey) + "\t" +
                                                std::to_string(besttHitToExtend.alnLength) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(besttHitToExtend.dbLen) + "\t" +
                                                std::to_string(besttHitToExtend.seqId) + "\t" +
                                                std::to_string(besttHitToExtend.rySeqId) + "\t" +
                                                "0\t" +  // Assuming this is a constant value
                                                std::to_string(query.length()) + "\t\n";

                        #pragma omp critical
                        {
                            if (outputFileTwo.is_open()) {
                                outputFileTwo << outputLine;
                            }
                        }
#endif

                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                    }

                }

                if (leftQueryOffset > 0 || rightQueryOffset > 0)
                  queryCouldBeExtended = true;

                if (!alnQueueReads.empty())
                    break;

                querySeqLen = query.length();
                querySeq = (char *) query.c_str();

                // update alignments

                for(size_t alnIdx = 0; alnIdx < tmpAlignmentsReads.size(); alnIdx++) {

                    unsigned int tId = sequenceDbr->getId(tmpAlignmentsReads[alnIdx].dbKey);
                    unsigned int tSeqLen = sequenceDbr->getSeqLen(tId);
                    char *tSeq = sequenceDbr->getData(tId, thread_idx);
                    bool deleteTargetSeq = false;
                    if (useReverse[tId])
                    {
                        tSeq = getNuclRevFragment(tSeq, tSeqLen, (NucleotideMatrix *) subMat);
                        deleteTargetSeq = true;
                    }
                        
                    int qStartPos = tmpAlignmentsReads[alnIdx].qStartPos;
                    int dbStartPos = tmpAlignmentsReads[alnIdx].dbStartPos;
                    int diag = (qStartPos + leftQueryOffset) - dbStartPos;

                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::ungappedAlignmentByDiagonal(
                                                                   querySeq, querySeqLen, tSeq, tSeqLen,
                                                                   diag, fastMatrix.matrix, par.rescoreMode);

                    updateNuclAlignment(tmpAlignmentsReads[alnIdx], alignment, querySeq, querySeqLen, tSeq, tSeqLen);

/*                     float ryId = rySeqId(tmpAlignments[alnIdx], querySeq, tSeq, nucleotideMap);
                    //std::cerr << "ryId " << ryId << std::endl;
                    tmpAlignments[alnIdx].rySeqId = ryId; */

                    if (deleteTargetSeq) {
                        delete[] tSeq;
                    }
                }


                // std::vector<unsigned int> maxAlign = getMaxAlnLen(tmpAlignmentsReads, queryKey);
                // maxAlnLeft = maxAlign[0];
                // maxAlnRight = maxAlign[1];
                // maxAlign.clear();

                consensus = consensusCaller(tmpAlignmentsReads, sequenceDbr, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat);
                updateSeqIdConsensusReads(tmpAlignmentsReads, sequenceDbr, consensus, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat, maxAlnLeft, maxAlnRight); 

                for(size_t alnIdx = 0; alnIdx < tmpAlignmentsReads.size(); alnIdx++) {
                    bool notInside = tmpAlignmentsReads[alnIdx].dbLen != tmpAlignmentsReads[alnIdx].alnLength;
                    const bool rightStart = tmpAlignmentsReads[alnIdx].dbStartPos == 0;
                    const bool leftStart = tmpAlignmentsReads[alnIdx].qStartPos == 0;
                    const bool isNotIdentity = (tmpAlignmentsReads[alnIdx].dbKey != queryKey);

                    if(tmpAlignmentsReads[alnIdx].seqId >= par.seqIdThr && (rightStart || leftStart) && isNotIdentity && notInside)
                        {
                            unsigned int tId = sequenceDbr->getId(tmpAlignmentsReads[alnIdx].dbKey);
                            char *tSeq = sequenceDbr->getData(tId, thread_idx);
                            //std::cerr << "b" << std::endl;
                            bool deleteTargetSeq = false;

                            if (tmpAlignmentsReads[alnIdx].isRevToAlignment){
                                tSeq = getNuclRevFragment(tSeq, tmpAlignmentsReads[alnIdx].dbLen, (NucleotideMatrix *) subMat);
                                deleteTargetSeq = true;
                            }

                            // #pragma omp critical
                            // std::cerr << "Before SECOND" << std::endl;
                            scorePerRes toAdd = r_s_pair(tmpAlignmentsReads[alnIdx], consensus, tSeq, querySeqLen,subDeamDiNuc, subDeamDiNucRev, maxAlnLeft, maxAlnRight, randAlnPenal, seqErrMatch, par.excessPenal);
                            // #pragma omp critical
                            // std::cerr << "After SECOND" << std::endl;
                            if ( toAdd.sRatio > par.likelihoodThreshold ) {
                                alnQueueReads.push( toAdd );
                            }
                            if (deleteTargetSeq) {
                                delete[] tSeq;
                            }
                        }
                }
            }

            if (queryCouldBeExtended)  {
                query.push_back('\n');
                __sync_or_and_fetch(&wasExtended[id], static_cast<unsigned char>(0x20));
                
                bool queryWasExtended = true;
                resultWriter.writeData(query.c_str(), query.size(), queryKey, thread_idx, queryWasExtended);
            }

            //contig.clear();
            notContig.clear();
        }


} // end parallel

// std::cerr << "TotalR" << "\t" << totalRight << std::endl;
// std::cerr << "countR" << "\t" << rightExtLongest << std::endl;
// std::cerr << "TotalL" << "\t" << totalLeft << std::endl;
// std::cerr << "countL" << "\t" << leftExtLongest << std::endl;
// std::cerr << "NumWasReverse\t" << wasReverse << "\tof\t" << totalQuer << std::endl; 

// add sequences that are not yet assembled
#pragma omp parallel for schedule(dynamic, 10000)
    for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        //bool couldExtend =  (wasExtended[id] & 0x10);
        bool isNotContig =  !(wasExtended[id] & 0x20);
        //bool wasNotUsed =  !(wasExtended[id] & 0x40);
        //bool wasNotExtended =  !(wasExtended[id] & 0x80);
        //bool wasUsed    =  (wasExtended[id] & 0x40);
        //if(isNotContig && wasNotExtended ){
        bool queryWasExtendedBefore = sequenceDbr->getExtData(id);
        if (isNotContig){
            char *querySeqData = sequenceDbr->getData(id, thread_idx);
            resultWriter.writeData(querySeqData, sequenceDbr->getEntryLen(id)-1, sequenceDbr->getDbKey(id), thread_idx, queryWasExtendedBefore);
        }
    }

    // cleanup
    resultWriter.close(true);
    alnReader->close();
    delete [] wasExtended;
    delete alnReader;
    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    sequenceDbr->close();
    delete sequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}


int ancientReadsResults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute assembly.\n";
    return doNuclAssembly1(par);
}

