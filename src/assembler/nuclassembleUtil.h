#include "LocalParameters.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include <libgab/libgab.h>


#include <limits>
#include <cstdint>
#include <queue>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <unordered_map>


#ifdef OPENMP
#include <omp.h>
#endif

typedef struct {
    int count[4][11];
} countDeamCov;

typedef struct {
    long double s[16];
} probSubstition;

typedef struct {
    long double p[4][4];
} diNucleotideProb;

typedef struct {
    long double s[12];
 } substitutionRates;

 struct scorePerRes {
    Matcher::result_t r;
    double sLenNorm = 0;
    double sRatio = 0; 
};

void updateNuclAlignment(Matcher::result_t &tmpAlignment, DistanceCalculator::LocalAlignment &alignment,
                                const char *querySeq, size_t querySeqLen, const char *tSeq, size_t tSeqLen);

void getSeqErrorProf(diNucleotideProb & seqErrMatch, long double err);

static void readNucSubstitionRatesFreq(const string filename,vector<substitutionRates> & subVec){
    igzstream subFP;

    subFP.open(filename.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (subFP.good()){
	vector<string> fields;
	string line;

	//header
	if ( !getline (subFP,line)){
	    std::cerr << "Unable to open file " << "\n";
        exit(1);
	}
	fields = allTokens(line,'\t');

	if(fields.size() != 12){
	    std::cerr << "Profile not 12 fields uniq1" << "\n";
        exit(1);
	}


	//probs
	while ( getline (subFP,line)){

	    fields = allTokens(line,'\t');

	    if(fields.size() != 12){
	    std::cerr << "Profile not 12 fields uniq2" << "\n";
        exit(1);
	    }

	    substitutionRates tempFreq;


	    for(unsigned int k=0;k<12;k++){
		//for(unsigned int t=0;t<=2;t++){
		tempFreq.s[k]=destringify<long double>(fields[k]);
		//}
	    }

	    subVec.emplace_back( tempFreq );
	}
	subFP.close();
    }else{
	    std::cerr << "Profile not 12 fields uniq3" << "\n";
        exit(1);
    }
}

//void createNoDamageMatrix(const std::string& filename);

void initDeamProbabilities(const std::string & deam5pfreqE,const std::string & deam3pfreqE, std::vector<substitutionRates> & sub5p, std::vector<substitutionRates> & sub3p, std::vector<diNucleotideProb> & allDeam, std::vector<diNucleotideProb> & revAllDeam);

std::vector<unsigned int> getMaxAlnLen(std::vector<Matcher::result_t> &alignments, unsigned int & queryKey);

class CompareNuclResultByScoreReads {
public:
    bool operator()(scorePerRes r1, scorePerRes r2) {
        double score1 = r1.sLenNorm;
        double score2 = r2.sLenNorm;

        if (score1 < score2) {
            return true;
        }
        return false;
    }
};

char* getNuclRevFragment(char* fragment, size_t fragLen, NucleotideMatrix *nuclMatrix);

float getRYSeqId(Matcher::result_t & res, char* querySeq,  char* targetSeq, std::unordered_map<char, int> & mapACGT);

void calcLikelihood(scorePerRes & scoredRes, char* querySeq, const char* targetSeq, std::vector<diNucleotideProb> & subDeamDiNuc, unsigned int maxAln, float randAlnPenal, diNucleotideProb & seq_err_match, float excessPenal);

std::pair<bool, bool> mgeFinder(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat, const std::vector<diNucleotideProb> & subDeamDiNuc, const std::vector<diNucleotideProb> & subDeamDiNucRev, const diNucleotideProb & seqErrMatch);

std::pair<bool, bool> mgeFinderContigs(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat);

std::string consensusCaller(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat);

void updateSeqIdConsensus(std::vector<Matcher::result_t> & alignments, DBReader<unsigned int>* sequenceDbr, std::string & consensus, char *querySeq, const unsigned int querySeqLen, unsigned int queryKey, unsigned int thread_idx, LocalParameters &par, NucleotideMatrix *subMat);

//bool calcLikelihoodCorrection(Matcher::result_t rightLongestCandi, Matcher::result_t candidate, std::string querySeq, std::string targetSeq, std::vector<diNucleotideProb> & subDeamDiNuc, diNucleotideProb & seqErrMatch, bool extSide);
std::vector<long double> calcLikelihoodCorrection(const Matcher::result_t & rightLongestCandi, const Matcher::result_t & candidate, const std::string & querySeq, const std::string & targetSeq, const std::vector<diNucleotideProb> & subDeamDiNuc, const diNucleotideProb & seqErrMatch, bool extSide);

int doNuclAssembly1(LocalParameters &par);

int doNuclAssembly2(LocalParameters &par);

// TODO:
// Function that calcualates the ancient specific mismatch count per result_t according to Johannes new method. 
// Inputs: result_t, damage matrix, sequence information

double deamMatches(Matcher::result_t res, int qBase, int tBase, unsigned int pos, unsigned int scoreAln, double matchLik);

float ancientMatchCount(Matcher::result_t res, char* querySeq, char* targetSeq, std::vector<diNucleotideProb> &subDeamDiNuc, std::unordered_map<char, int> &nucleotideMap);

float ancientMatchCountContig(Matcher::result_t res, char* querySeq, char* targetSeq, std::vector<diNucleotideProb> &subDeamDiNuc, std::unordered_map<char, int> &nucleotideMap);

float getSubSeqId(Matcher::result_t &aln, const char *querySeq, const char *tSeq);
 

