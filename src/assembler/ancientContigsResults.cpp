#include "LocalParameters.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include "nuclassembleUtil.h"

#include <limits>
#include <cstdint>
#include <queue>
#include <vector>

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

//#define DEBUGCOV


class CompareNuclResultByScoreContigs {
public:
    // bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
        
    //     float mm_count1 = r1.alnLength - r1.deamMatch;
    //     float mm_count2 = r2.alnLength - r2.deamMatch;

    //     float alpha1 = mm_count1 + 1;
    //     float alpha2 = mm_count2 + 1;
    //     float beta1 = r1.deamMatch + 1;
    //     float beta2 = r2.deamMatch + 1;

    //     //double c=(std::tgamma(beta1+beta2)*std::tgamma(alpha1+beta2))/(std::tgamma(alpha1+beta1+beta2)*std::tgamma(beta1));
    //     double log_c = (std::lgamma(beta1+beta2)+std::lgamma(alpha1+beta1))-(std::lgamma(alpha1+beta1+beta2)+std::lgamma(beta1));

    //     //double r = 1.0; // r_0 =1
    //     double log_r = 0.0;
    //     double p = 0.0;
    //     for (size_t idx = 0; idx < alpha2; idx++) {

    //         p += exp(log_r + log_c);
    //         //r *= ((alpha1+idx)*(beta2+idx))/((idx+1)*(idx+alpha1+beta1+beta2));
    //         log_r = log(alpha1+idx)+log(beta2+idx)-(log(idx+1) + log(idx+alpha1+beta1+beta2)) + log_r;
    //     }

    bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
        unsigned int mm_count1 = (1 - r1.seqId) * r1.alnLength + 0.5;
        unsigned int mm_count2 = (1 - r2.seqId) * r2.alnLength + 0.5;

        unsigned int alpha1 = mm_count1 + 1;
        unsigned int alpha2 = mm_count2 + 1;
        unsigned int beta1 = r1.alnLength - mm_count1 + 1;
        unsigned int beta2 = r2.alnLength - mm_count2 + 1;

        //double c=(std::tgamma(beta1+beta2)*std::tgamma(alpha1+beta2))/(std::tgamma(alpha1+beta1+beta2)*std::tgamma(beta1));
        double log_c = (std::lgamma(beta1+beta2)+std::lgamma(alpha1+beta1))-(std::lgamma(alpha1+beta1+beta2)+std::lgamma(beta1));

        //double r = 1.0; // r_0 =1
        double log_r = 0.0;
        double p = 0.0;
        for (size_t idx = 0; idx < alpha2; idx++) {

            p += exp(log_r + log_c);
            //r *= ((alpha1+idx)*(beta2+idx))/((idx+1)*(idx+alpha1+beta1+beta2));
            log_r = log(alpha1+idx)+log(beta2+idx)-(log(idx+1) + log(idx+alpha1+beta1+beta2)) + log_r;
        }
        //p *= c;

        if (p < 0.45)
            return true;
        if (p  > 0.55)
            return false;
        // if ( r1.alnLength > 200 && r2.alnLength > 200 ){
        //     if ( r1.seqId > r2.seqId ){
        //         return false;
        //     }
        //     if ( r2.seqId > r1.seqId ){
        //         return true;
        //     }
        // }
        if (r1.alnLength < r2.alnLength)
            return true;
        if (r1.alnLength > r2.alnLength)
            return false;
        // before simply longer alignment


        // if ( r1.dbLen < 200 || r2.dbLen < 200 ){
        //     if (r1.alnLength < r2.alnLength)
        //         return true;
        //     if (r1.alnLength > r2.alnLength)
        //         return false;
        //     if ( r1.seqId > r2.seqId ){
        //         return false;
        //     }
        //     if ( r2.seqId > r1.seqId ){
        //         return true;
        //     }
        // }
        // if ( r1.seqId > r2.seqId ){
        //     return false;
        // }
        // if ( r2.seqId > r1.seqId ){
        //     return true;
        // }
        // if (r1.alnLength > r2.alnLength)
        //     return false;
        // if (r1.alnLength < r2.alnLength)
        //     return true;



        // if (r1.dbLen < r2.dbLen)
        //     return true;
        // if (r1.dbLen > r2.dbLen)
        //     return false;

        return true;
    }
};

typedef std::priority_queue<Matcher::result_t, std::vector<Matcher::result_t> , CompareNuclResultByScoreContigs> QueueByScoreNuclContigs;
Matcher::result_t selectNuclFragmentToExtendContigs(QueueByScoreNuclContigs &alignments,
                                             unsigned int queryKey) {
    // results are ordered by score
    while (alignments.empty() == false){
        Matcher::result_t res = alignments.top();
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);

        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}



int doNuclAssembly2(LocalParameters &par) {
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
    std::cerr << userInput << std::endl;
    
    std::string high5 = userInput + "5p.prof";
    std::string high3 = userInput + "3p.prof";

#ifdef DEBUGCOV
    auto currentTime = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(currentTime);

    std::stringstream one;
    one << "all_candidates_" << time << ".tsv";
    std::string filenameOne = one.str();
    std::ofstream outputFileAll;  // Declare the ofstream object globally if possible
    outputFileAll.open(filenameOne, std::ios::app);  // Open once, append mode

    std::stringstream two;
    two << "ext_candidates_" << time << ".tsv";
    std::string filenameTwo = two.str();
    std::ofstream outputFileExt;  // Declare the ofstream object globally if possible
    outputFileExt.open(filenameTwo, std::ios::app);  // Open once, append mode

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

        // Vector for substitution rates given that we alwys use first 5, last 5 and mid position
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

// #ifdef DEBUGCOV
//         // debugging coverage
//         //std::vector<std::vector<double>> allCovsDebug;
//         // std::vector<std::vector<double>> besthitCoverages; 
// #endif


#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            progress.updateProgress();

// #ifdef DEBUGCOV
//             //debugging coverage
//             std::vector<double> dbCovs;
// #endif

            unsigned int queryKey = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id, thread_idx);
            unsigned int querySeqLen = sequenceDbr->getSeqLen(id);
            std::string query(querySeq, querySeqLen); // no /n/0

            char *alnData = alnReader->getDataByDBKey(queryKey, thread_idx);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);

            bool queryCouldBeExtended = false;
            QueueByScoreNuclContigs alnQueue;

            bool mgeFoundRight = false;
            bool mgeFoundLeft = false;

            //counter for left and right extensions
            unsigned int countRightExt = 0;
            unsigned int countLeftExt = 0;

            std::vector<Matcher::result_t> contigs;
            contigs.reserve(300);

            // get the aligned positions right so it stays the same for all following procedures and does not to be done anymore:




            for (unsigned int idx = 0; idx < alignments.size(); idx++){
                // now retrieve the other candidate extensions and get their sequences

                unsigned int aln2updateId = sequenceDbr->getId(alignments[idx].dbKey);
                unsigned int aln2updateLen = sequenceDbr->getSeqLen(aln2updateId);

                char *aln2updateSequ = sequenceDbr->getData(aln2updateId, thread_idx);
                std::string aln2updateSeq;

                int rawScore = static_cast<int>(evaluer.computeRawScoreFromBitScore(alignments[idx].score) + 0.5);
                float scorePerCol = static_cast<float>(rawScore) / static_cast<float>(alignments[idx].alnLength + 0.5);

                if (alignments[idx].qStartPos > alignments[idx].qEndPos) {
                    // Convert the reversed fragment to std::string
                    useReverse[sequenceDbr->getId(alignments[idx].dbKey)] = true;
                    alignments[idx].isRevToAlignment = true;

                    std::swap(alignments[idx].qStartPos, alignments[idx].qEndPos);
                    unsigned int dbStartPos = alignments[idx].dbStartPos;
                    alignments[idx].dbStartPos = alignments[idx].dbLen - alignments[idx].dbEndPos - 1;
                    alignments[idx].dbEndPos= alignments[idx].dbLen - dbStartPos - 1;

                    char *rightExtCandiSeqTmp = getNuclRevFragment(aln2updateSequ, aln2updateLen, (NucleotideMatrix *)subMat);
                    aln2updateSeq = std::string(rightExtCandiSeqTmp, aln2updateLen);
                    delete[] rightExtCandiSeqTmp;
                } else {
                    // Directly convert the raw sequence to std::string
                    aln2updateSeq = std::string(aln2updateSequ, aln2updateLen);
                    useReverse[sequenceDbr->getId(alignments[idx].dbKey)] = false;
                    alignments[idx].isRevToAlignment = false;
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
                alignments[idx].seqId = seqId;
                alignments[idx].rySeqId = rySeqId;

                if ( alignments[idx].seqId >= par.mergeSeqIdThr && alignments[idx].rySeqId >= par.rySeqIdThr){
                    contigs.push_back(alignments[idx]); 
                }
            }


            if ( contigs.size() > 1){
                std::string consensus = consensusCaller(contigs, sequenceDbr, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat);
                updateSeqIdConsensus(contigs, sequenceDbr, consensus, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat); 
            }
            alignments.clear();

            // std::pair<bool, bool> mgefound = mgeFinderContigs(contigs, sequenceDbr, querySeq, querySeqLen, queryKey, thread_idx, par, (NucleotideMatrix *) subMat);
            // mgeFoundLeft = mgefound.first;
            // mgeFoundRight = mgefound.second;


            // fill queue
            for (size_t alnIdx = 0; alnIdx < contigs.size(); alnIdx++) {

                if ( contigs[alnIdx].seqId >= par.mergeSeqIdThr && contigs[alnIdx].rySeqId >= par.rySeqIdThr )
                {
                    // std::vector<diNucleotideProb> subDeamDiNucRef = contigs[alnIdx].isRevToAlignment ? subDeamDiNucRev : subDeamDiNuc;
                    // float ancientMatches = ancientMatchCount(contigs[alnIdx], querySeq, tSeq, subDeamDiNucRef, nucleotideMap);
                    // contigs[alnIdx].deamMatch = ancientMatches;

                    alnQueue.push(contigs[alnIdx]);
                    if (contigs.size() > 1){
                        __sync_or_and_fetch(&wasExtended[sequenceDbr->getId(contigs[alnIdx].dbKey)],
                                            static_cast<unsigned char>(0x40));
                    }
                }
            }



            // Extension process begins here
            std::vector<Matcher::result_t> tmpAlignments;
            tmpAlignments.reserve(contigs.size());

            while (!alnQueue.empty()) {

                unsigned int leftQueryOffset = 0;
                unsigned int rightQueryOffset = 0;
                tmpAlignments.clear();
                Matcher::result_t besttHitToExtend;

                if ( mgeFoundLeft && mgeFoundRight ){
                    QueueByScoreNuclContigs alnQueue;
                    break;
                } 


                while ((besttHitToExtend = selectNuclFragmentToExtendContigs(alnQueue, queryKey)).dbKey != UINT_MAX) {

                    unsigned int targetId = sequenceDbr->getId(besttHitToExtend.dbKey);
                    if (targetId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Could not find targetId  " << besttHitToExtend.dbKey
                                            << " in database " << sequenceDbr->getDataFileName() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *targetSeq = sequenceDbr->getData(targetId, thread_idx);
                    unsigned int targetSeqLen = sequenceDbr->getSeqLen(targetId);

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


#ifdef DEBUGCOV
                    std::string targetPrintAll(targetSeq, targetSeqLen); 
                    std::string outputLineAll = std::to_string(queryKey) + "\t" +
                                                std::to_string(besttHitToExtend.alnLength) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(besttHitToExtend.dbLen) + "\t" +
                                                std::to_string(besttHitToExtend.seqId) + "\t" +
                                                std::to_string(besttHitToExtend.rySeqId) + "\t" +
                                                "2\t" +  // Assuming this is a constant value
                                                std::to_string(query.length()) + "\t" +
                                                std::to_string(besttHitToExtend.qStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.qEndPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbEndPos) + "\t" + 
                                                query + "\t" +
                                                targetPrintAll + "\t\n";

                        #pragma omp critical
                        {
                            if (outputFileAll.is_open()) {
                                outputFileAll << outputLineAll;
                            }
                        }
#endif

                    bool solidRightExt = !(countRightExt == 1 && besttHitToExtend.alnLength <= 100);
                    bool solidLeftExt = !(countLeftExt == 1 && besttHitToExtend.alnLength <= 100);

                    if (dbStartPos == 0 && qEndPos == (querySeqLen - 1) && solidRightExt && !mgeFoundRight) {
                        //right extension

                        if (rightQueryOffset > 0) {
                            tmpAlignments.push_back(besttHitToExtend);
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

#ifdef DEBUGCOV
                        bool wasReverse = useReverse[targetId];
                        std::string queryOld = query;
                        std::string result = query + fragment;
#endif

                        query += fragment;
                        rightQueryOffset += fragLen;

#ifdef DEBUGCOV
                    #pragma omp critical
                    {
                    std::string targetPrint(targetSeq, targetSeqLen); 
                    std::string outputLine = std::to_string(queryKey) + "\t" +
                                                std::to_string(besttHitToExtend.alnLength) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(besttHitToExtend.dbLen) + "\t" +
                                                std::to_string(besttHitToExtend.seqId) + "\t" +
                                                std::to_string(besttHitToExtend.rySeqId) + "\t" +
                                                "1\t" +  // Assuming this is a constant value
                                                std::to_string(query.length()) + "\t" +
                                                std::to_string(besttHitToExtend.qStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.qEndPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbEndPos) + "\t" + 
                                                std::to_string(wasReverse) + "\t" + 
                                                queryOld + "\t" +
                                                targetPrint + "\t" +
                                                result + "\t\n";

                        if (outputFileExt.is_open()) {
                            outputFileExt << outputLine;
                        }
                    }
#endif

                    //update that dbKey was used in assembly
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                    }
                    else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1) && solidLeftExt && !mgeFoundLeft) {
                        //left extension

                        if(leftQueryOffset > 0) {
                            tmpAlignments.push_back(besttHitToExtend);
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

#ifdef DEBUGCOV
                        bool wasReverse = useReverse[targetId];
                        std::string queryOld = query;
                        std::string result = query + fragment;
#endif

                        query = fragment + query;
                        leftQueryOffset += fragLen;

#ifdef DEBUGCOV
                    #pragma omp critical
                    {
                    std::string targetPrint(targetSeq, targetSeqLen); 
                    std::string outputLine = std::to_string(queryKey) + "\t" +
                                                std::to_string(besttHitToExtend.alnLength) + "\t" +
                                                std::to_string(querySeqLen) + "\t" +
                                                std::to_string(besttHitToExtend.dbLen) + "\t" +
                                                std::to_string(besttHitToExtend.seqId) + "\t" +
                                                std::to_string(besttHitToExtend.rySeqId) + "\t" +
                                                "0\t" +  // Assuming this is a constant value
                                                std::to_string(query.length()) + "\t" +
                                                std::to_string(besttHitToExtend.qStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.qEndPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbStartPos) + "\t" + 
                                                std::to_string(besttHitToExtend.dbEndPos) + "\t" + 
                                                std::to_string(wasReverse) + "\t" + 
                                                queryOld + "\t" +
                                                targetPrint + "\t" +
                                                result + "\t\n";


                        if (outputFileExt.is_open()) {
                            outputFileExt << outputLine;
                        }
                    }
#endif
                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                    }
                }

                if (leftQueryOffset > 0 || rightQueryOffset > 0)
                  queryCouldBeExtended = true;

                if (!alnQueue.empty())
                    break;

                querySeqLen = query.length();
                querySeq = (char *) query.c_str();
    
                // Re-initilaize the extension count
                countRightExt = 0;
                countLeftExt = 0;

                // update alignments
                for(size_t alnIdx = 0; alnIdx < tmpAlignments.size(); alnIdx++) {

                    unsigned int tId = sequenceDbr->getId(tmpAlignments[alnIdx].dbKey);
                    unsigned int tSeqLen = sequenceDbr->getSeqLen(tId);
                    char *tSeq = sequenceDbr->getData(tId, thread_idx);
                    if (useReverse[tId])
                        tSeq = getNuclRevFragment(tSeq, tSeqLen, (NucleotideMatrix *) subMat);

                    int qStartPos = tmpAlignments[alnIdx].qStartPos;
                    int dbStartPos = tmpAlignments[alnIdx].dbStartPos;
                    int diag = (qStartPos + leftQueryOffset) - dbStartPos;

                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::ungappedAlignmentByDiagonal(
                                                                   querySeq, querySeqLen, tSeq, tSeqLen,
                                                                   diag, fastMatrix.matrix, par.rescoreMode);

                    updateNuclAlignment(tmpAlignments[alnIdx], alignment, querySeq, querySeqLen, tSeq, tSeqLen);

                    float ryId = getRYSeqId(tmpAlignments[alnIdx], querySeq,  tSeq, ryMap);
                    tmpAlignments[alnIdx].rySeqId = ryId;

                    // refill queue
                    //if(tmpAlignments[alnIdx].seqId >= par.mergeSeqIdThr && tmpAlignments[alnIdx].rySeqId >= par.rySeqIdThr){
                    if(tmpAlignments[alnIdx].seqId >= par.mergeSeqIdThr && tmpAlignments[alnIdx].rySeqId >= par.rySeqIdThr){ 
                        std::vector<diNucleotideProb> subDeamDiNucRef = tmpAlignments[alnIdx].isRevToAlignment ? subDeamDiNucRev : subDeamDiNuc;
                        float ancientMatches = ancientMatchCount(tmpAlignments[alnIdx], querySeq, tSeq, subDeamDiNucRef, nucleotideMap);
                        tmpAlignments[alnIdx].deamMatch = ancientMatches;
                        alnQueue.push(tmpAlignments[alnIdx]);

                        size_t tmpDbKey = tmpAlignments[alnIdx].dbKey;
                        const bool notRightStartAndLeftStart = !(tmpAlignments[alnIdx].dbStartPos == 0 && tmpAlignments[alnIdx].qStartPos == 0 );
                        const bool rightStart = tmpAlignments[alnIdx].dbStartPos == 0 && (tmpAlignments[alnIdx].dbEndPos != static_cast<int>(tmpAlignments[alnIdx].dbLen)-1);
                        const bool leftStart = tmpAlignments[alnIdx].qStartPos == 0 && (tmpAlignments[alnIdx].qEndPos != static_cast<int>(tmpAlignments[alnIdx].qLen)-1);
                        const bool isNotIdentity = (tmpDbKey != queryKey);

                        // distinguish between right and left here
                        unsigned int dbStartPos = tmpAlignments[alnIdx].dbStartPos;
                        unsigned int dbEndPos = tmpAlignments[alnIdx].dbEndPos;
                        unsigned int qStartPos = tmpAlignments[alnIdx].qStartPos;
                        unsigned int qEndPos = tmpAlignments[alnIdx].qEndPos;               

                        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){

                            // count ALL right and left extensions; needed for pseudo stable kmer filter
                            if (dbStartPos == 0 && qEndPos == (querySeqLen - 1)) {
                                // right extension
                                countRightExt++; 
                            }
                            else if (qStartPos == 0 && dbEndPos == (tmpAlignments[alnIdx].dbLen - 1)) {
                                // left extension
                                countLeftExt++;
                            }
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

            contigs.clear();
        }

    } // end parallel

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

int ancientContigsResults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute assembly.\n";
    return doNuclAssembly2(par);
}