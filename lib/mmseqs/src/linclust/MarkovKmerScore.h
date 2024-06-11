//
// Created by mad on 2019-02-27.
//

#ifndef MMSEQS_KMERMARKOVSCORE_H
#define MMSEQS_KMERMARKOVSCORE_H

#include <Indexer.h>

namespace MarkovScores {
    extern const int MARKOV_ORDER;
    extern const float MEDIAN_SCORE;
    extern float markov5Scores[];
}
class MarkovKmerScore{
public:



    static float scoreKmer(const unsigned char * kmer, unsigned char kmerSize){
        float totalSocore = 0.0;
        for(int pos = 0; pos < kmerSize - MarkovScores::MARKOV_ORDER; pos++){
            size_t lookupIdx = Indexer::computeKmerIdx(&kmer[pos], MarkovScores::MARKOV_ORDER+1);
            totalSocore += MarkovScores::markov5Scores[lookupIdx];
        }
        return totalSocore;
    }

    static int adjustedLength(const unsigned char * kmer, unsigned char kmerSize, float minScoreThr){
        float totalSocore = 0.0;
        int pos = 0;
        while(totalSocore < minScoreThr && pos < kmerSize - MarkovScores::MARKOV_ORDER){
            size_t lookupIdx = Indexer::computeKmerIdx(&kmer[pos], MarkovScores::MARKOV_ORDER+1);
            float score = MarkovScores::markov5Scores[lookupIdx];
            totalSocore += score;
            pos++;
        }
        return pos+MarkovScores::MARKOV_ORDER;
    }
};

#endif //MMSEQS_KMERMARKOVSCORE_H
