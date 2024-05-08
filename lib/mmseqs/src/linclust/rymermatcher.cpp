// include xxhash early to avoid incompatibilites with SIMDe
#define XXH_INLINE_ALL
#include "xxhash.h"

#include "kmermatcher.h"
#include "Debug.h"
#include "Indexer.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "NucleotideMatrix.h"
#include "tantan.h"
#include "QueryMatcher.h"
#include "KmerGenerator.h"
#include "MarkovKmerScore.h"
#include "FileUtil.h"
#include "FastSort.h"

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <limits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif
#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

uint64_t hashUInt64(uint64_t in, uint64_t seed, bool rymatch) {
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    in = __builtin_bswap64(in);
#endif
    return XXH64(&in, sizeof(uint64_t), seed);
}

template <int TYPE, typename T>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution, bool rymatch){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.rymerSize;
    ProbabilityMatrix *probMatrix = NULL;
    if (par.maskMode == 1) {
        probMatrix = new ProbabilityMatrix(*subMat);
    }

    ScoreMatrix two;
    ScoreMatrix three;
    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        two = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(seqDbr.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned short * scoreDist= new unsigned short[65536];
        unsigned int * hierarchicalScoreDist= new unsigned int[128];

        const int adjustedKmerSize = (par.adjustKmerLength) ? std::min( par.rymerSize+5, 23) :   par.rymerSize;
        Sequence seq(par.maxSeqLen, querySeqType, subMat, adjustedKmerSize, par.spacedKmer, false, true, par.spacedKmerPattern);
        KmerGenerator* generator;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            generator = new KmerGenerator( par.rymerSize, subMat->alphabetSize, 150);
            generator->setDivideStrategy(&three, &two);
        }
        Indexer idxer(subMat->alphabetSize - 1,  par.rymerSize);
        const unsigned int BUFFER_SIZE = 1048576;
        size_t bufferPos = 0;
        KmerPosition<T> * threadKmerBuffer = new KmerPosition<T>[BUFFER_SIZE];
        SequencePosition * kmers = (SequencePosition *) malloc((par.pickNbest * (par.maxSeqLen + 1) + 1) * sizeof(SequencePosition));
        size_t kmersArraySize = par.maxSeqLen;
        const size_t flushSize = 100000000;
        size_t iterations = static_cast<size_t>(ceil(static_cast<double>(seqDbr.getSize()) / static_cast<double>(flushSize)));
        for (size_t i = 0; i < iterations; i++) {
            size_t start = (i * flushSize);
            size_t bucketSize = std::min(seqDbr.getSize() - (i * flushSize), flushSize);

#pragma omp for schedule(dynamic, 100)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                memset(scoreDist, 0, sizeof(unsigned short) * 65536);
                memset(hierarchicalScoreDist, 0, sizeof(unsigned int) * 128);

                seq.mapSequence(id, seqDbr.getDbKey(id), seqDbr.getData(id, thread_idx), seqDbr.getSeqLen(id));

                size_t seqHash =  SIZE_T_MAX;
                //TODO, how to handle this in reverse?
                if(hashWholeSequence){
                    seqHash = Util::hash(seq.numSequence, seq.L);
                    seqHash = hashUInt64(seqHash, par.hashShift, true);
                }

                maskSequence(par.maskMode, par.maskLowerCaseMode, seq, subMat->aa2num[static_cast<int>('X')], probMatrix);

                size_t seqKmerCount = 0;
                unsigned int seqId = seq.getDbKey();
                while (seq.hasNextKmer()) {
                    unsigned char *kmer = (unsigned char*) seq.nextKmer();
                    if(seq.kmerContainsX()){
                        continue;
                    }
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        NucleotideMatrix * nuclMatrix = (NucleotideMatrix*)subMat;
                        size_t kmerLen =  par.rymerSize;

                        // lOUIS WAS HERE
                        //if (par.krymatcher.aminoacids)
                        size_t kmerIdx = Indexer::computeRYmerIdx(kmer, kmerLen);
                        size_t revkmerIdx = Util::revComplement(kmerIdx, kmerLen);
                        // skip forward and rev. identical k-mers.
                        // We can not know how to align these afterwards
                        if(revkmerIdx == kmerIdx){
                            continue;
                        }
                        bool pickReverseKmer = (revkmerIdx<kmerIdx);
                        kmerIdx = (pickReverseKmer) ? revkmerIdx : kmerIdx;
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift, true);

                        if(par.adjustKmerLength) {
                            unsigned char revKmer[32];
                            unsigned char * kmerToHash = kmer;
                            if(pickReverseKmer){
                                for(int pos = static_cast<int>(adjustedKmerSize)-1; pos > -1; pos--){
                                    revKmer[(adjustedKmerSize - 1) - pos]=nuclMatrix->reverseResidue(kmer[pos]);
                                }
                                kmerToHash = revKmer;
                            }
                            kmerLen = MarkovKmerScore::adjustedLength(kmerToHash, adjustedKmerSize,
                                                                      (par.rymerSize - MarkovScores::MARKOV_ORDER) * MarkovScores::MEDIAN_SCORE);
                            longestKmer = std::max(kmerLen, longestKmer);
                            kmerIdx = Indexer::computeRYmerIdx(kmerToHash, kmerLen);
                        }

                        // set signed bit for normal kmers to make the  SIZE_T_MAX logic easier
                        // reversed kmers do not have a signed bit
                        size_t kmerRev = (pickReverseKmer) ? BIT_CLEAR(kmerIdx, 63) : BIT_SET(kmerIdx, 63);
                        (kmers + seqKmerCount)->kmer = kmerRev;
                        int pos = seq.getCurrentPosition();
                        (kmers + seqKmerCount)->pos = (pickReverseKmer) ? (seq.L) - pos - kmerLen : pos;
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
                        seqKmerCount++;
                    } else if(TYPE == Parameters::DBTYPE_HMM_PROFILE) {
                        std::pair<size_t*, size_t>  scoreMat = generator->generateKmerList(kmer, true);
//                        std::cout << scoreMat.elementSize << std::endl;
                        for(size_t kmerPos = 0; kmerPos < scoreMat.second && kmerPos < static_cast<size_t >(par.pickNbest); kmerPos++){
                            size_t kmerIdx = scoreMat.first[kmerPos];
                            (kmers + seqKmerCount)->kmer = kmerIdx;
                            (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                            const unsigned short hash = hashUInt64(kmerIdx, par.hashShift, true);
                            (kmers + seqKmerCount)->score = hash;
                            scoreDist[hash]++;
                            hierarchicalScoreDist[hash >> 9]++;
                            seqKmerCount++;
                        }
                    } else {
                        size_t kmerIdx = idxer.int2index(kmer, 0, par.rymerSize);
                        (kmers + seqKmerCount)->kmer = kmerIdx;
                        (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift, true);
//                        (kmers + seqKmerCount)->score = hash;
//                        const unsigned short hash = circ_hash(kmer, par.rymerSize, 5);
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
//                        std::cout << seqId << "\t" << (kmers + seqKmerCount)->score << "\t" << (kmers + seqKmerCount)->pos << std::endl;

                        seqKmerCount++;
                    }
                    if(seqKmerCount >= kmersArraySize){
                        kmersArraySize = seq.getMaxLen();
                        kmers = (SequencePosition *) realloc(kmers, (par.pickNbest * (kmersArraySize + 1) + 1) * sizeof(SequencePosition));
                    }

                }
                float kmersPerSequenceScale = (TYPE == Parameters::DBTYPE_NUCLEOTIDES) ? par.kmersPerSequenceScale.nucleotides
                                                                                       : par.kmersPerSequenceScale.aminoacids;
                size_t kmerConsidered = std::min(static_cast<size_t >(par.kmersPerSequence  - 1 + (kmersPerSequenceScale * seq.L)), seqKmerCount);

                unsigned int threshold = 0;
                size_t kmerInBins = 0;
                if (seqKmerCount > 0) {
                    size_t hierarchicaThreshold = 0;
                    for(hierarchicaThreshold = 0; hierarchicaThreshold < 128 && kmerInBins < kmerConsidered; hierarchicaThreshold++){
                        kmerInBins += hierarchicalScoreDist[hierarchicaThreshold];
                    }
                    hierarchicaThreshold -= (hierarchicaThreshold > 0) ? 1: 0;
                    kmerInBins -= hierarchicalScoreDist[hierarchicaThreshold];
                    for(threshold = hierarchicaThreshold*512; threshold <= USHRT_MAX && kmerInBins < kmerConsidered; threshold++){
                        kmerInBins += scoreDist[threshold];
                    }
                }
                int tooMuchElemInLastBin = (kmerInBins - kmerConsidered);

                // add k-mer to represent the identity
                if (static_cast<unsigned short>(seqHash) >= hashStartRange && static_cast<unsigned short>(seqHash) <= hashEndRange) {
                    threadKmerBuffer[bufferPos].kmer = seqHash;
                    threadKmerBuffer[bufferPos].id = seqId;
                    threadKmerBuffer[bufferPos].pos = 0;
                    threadKmerBuffer[bufferPos].seqLen = seq.L;
                    if(hashDistribution != NULL){
                        __sync_fetch_and_add(&hashDistribution[static_cast<unsigned short>(seqHash)], 1);
                    }
                    bufferPos++;
                    if (bufferPos >= BUFFER_SIZE) {
                        size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                        if(writeOffset + bufferPos < kmerArraySize){
                            if(kmerArray!=NULL){
                                memcpy(kmerArray + writeOffset, threadKmerBuffer, sizeof(KmerPosition<T>) * bufferPos);
                            }
                        } else{
                            Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset="<< writeOffset
                                                << ", kmerBufferPos=" << bufferPos
                                                << ", kmerArraySize=" << kmerArraySize <<".\n";
                            EXIT(EXIT_FAILURE);
                        }
                        bufferPos = 0;
                    }
                }

                if(par.ignoreMultiKmer){
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScoreReverse);
                    }else{
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScore);
                    }
                }
                size_t selectedKmer = 0;
                for (size_t kmerIdx = 0; kmerIdx < seqKmerCount && selectedKmer < kmerConsidered; kmerIdx++) {

                    /* skip repeated kmer */
                    if (par.ignoreMultiKmer) {
                        size_t kmer = (kmers + kmerIdx)->kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmerIdx + 1 < seqKmerCount) {
                            size_t nextKmer = (kmers + kmerIdx + 1)->kmer;
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                nextKmer = BIT_SET(nextKmer, 63);
                            }
                            if (kmer == nextKmer) {
                                while (kmer == nextKmer && kmerIdx < seqKmerCount) {
                                    kmerIdx++;
                                    if(kmerIdx >= seqKmerCount)
                                        break;
                                    nextKmer = (kmers + kmerIdx)->kmer;
                                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                        nextKmer = BIT_SET(nextKmer, 63);
                                    }
                                }
                            }
                        }
                        if(kmerIdx >= seqKmerCount)
                            break;
                    }

                    if ((kmers + kmerIdx)->score < threshold ){
                        // this if is needed to avoid extracting too much elements in the last bin
                        if((kmers + kmerIdx)->score == (threshold - 1) && tooMuchElemInLastBin){
                            tooMuchElemInLastBin--;
                            threshold -= (tooMuchElemInLastBin == 0) ? 1 : 0;
                        }
//                        std::cout << seqId << "\t" << (kmers + kmerIdx)->score << "\t" << (kmers + kmerIdx)->pos << std::endl;

                        selectedKmer++;
                        if ((kmers + kmerIdx)->score >= hashStartRange && (kmers + kmerIdx)->score <= hashEndRange)
                        {
//                            {
//                                size_t tmpKmerIdx= (kmers + kmerIdx)->kmer;
//                                tmpKmerIdx=BIT_CLEAR(tmpKmerIdx, 63);
//                                std::cout << seqId << "\t" << (kmers + kmerIdx)->score << "\t" << tmpKmerIdx << std::endl;
//                            }
                            threadKmerBuffer[bufferPos].kmer = (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].seqLen = seq.L;
                            bufferPos++;
                            if(hashDistribution != NULL){
                                __sync_fetch_and_add(&hashDistribution[(kmers + kmerIdx)->score], 1);
                            }

                            if (bufferPos >= BUFFER_SIZE) {
                                size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                                if(writeOffset + bufferPos < kmerArraySize){
                                    if(kmerArray!=NULL) {
                                        memcpy(kmerArray + writeOffset, threadKmerBuffer,
                                               sizeof(KmerPosition<T>) * bufferPos);
                                    }
                                } else{
                                    Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset="<< writeOffset
                                                        << ", kmerBufferPos=" << bufferPos
                                                        << ", kmerArraySize=" << kmerArraySize <<".\n";

                                    EXIT(EXIT_FAILURE);
                                }

                                bufferPos = 0;
                            }
                        }
                    }
                }
            }
#pragma omp barrier
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            if (thread_idx == 0) {
                seqDbr.remapData();
            }
#pragma omp barrier
        }


        if(bufferPos > 0){
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            if(kmerArray != NULL){
                memcpy(kmerArray+writeOffset, threadKmerBuffer, sizeof(KmerPosition<T>) * bufferPos);
            }
        }
        free(kmers);
        delete[] threadKmerBuffer;
        delete[] hierarchicalScoreDist;
        delete[] scoreDist;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            delete generator;
        }
    }

    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(three);
        ExtendedSubstitutionMatrix::freeScoreMatrix(two);
    }

    if (probMatrix != NULL) {
        delete probMatrix;
    }
    return std::make_pair(offset, longestKmer);
}

template <typename T>
KmerPosition<T> * doComputation(size_t totalKmers, size_t hashStartRange, size_t hashEndRange, std::string splitFile,
                                DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix  * subMat, bool rymatch) {

    KmerPosition<T> * hashSeqPair = initKmerPositionMemory<T>(totalKmers);
    size_t elementsToSort;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL, true);
        elementsToSort = ret.first;
        par.rymerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.rymerSize << "\n";
    }else{
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL, true);
        elementsToSort = ret.first;
    }
    if(hashEndRange == SIZE_T_MAX){
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T>::compareRepSequenceAndIdAndPosReverse);
    }else{
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    // assign rep. sequence to same kmer members
    // The longest sequence is the first since we sorted by kmer, seq.Len and id
    size_t writePos;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr);
    }else{
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr);
    }

    // sort by rep. sequence (stored in kmer) and sequence id
    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T>::compareRepSequenceAndIdAndDiagReverse);
    }else{
        SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T>::compareRepSequenceAndIdAndDiag);
    }
    //kx::radix_sort(hashSeqPair, hashSeqPair + elementsToSort, SequenceComparision());
//    for(size_t i = 0; i < writePos; i++){
//        std::cout << BIT_CLEAR(hashSeqPair[i].kmer, 63) << "\t" << hashSeqPair[i].id << "\t" << hashSeqPair[i].pos << std::endl;
//    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    if(hashEndRange != SIZE_T_MAX){
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T>(splitFile, hashSeqPair, writePos + 1);
        }else{
            writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T>(splitFile, hashSeqPair, writePos + 1);
        }
        delete [] hashSeqPair;
        hashSeqPair = NULL;
    }
    return hashSeqPair;
}


template <typename T>
int rymermatcherInner(Parameters& par, DBReader<unsigned int>& seqDbr) {

    int querySeqType = seqDbr.getDbtype();
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }else {
        if (par.alphabetSize.aminoacids == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.aminoacids, 8.0, -0.2f);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.aminoacids, 2.0);
        }
    }

    //seqDbr.readMmapedDataInMemory();

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.nucleotides : par.kmersPerSequenceScale.aminoacids;
    size_t totalKmers = computeKmerCount(seqDbr, par.rymerSize, par.kmersPerSequence, kmersPerSequenceScale);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T>(totalKmers);
    // compute splits
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
    size_t totalKmersPerSplit = std::max(static_cast<size_t>(1024+1),
                                         static_cast<size_t>(std::min(totalSizeNeeded, memoryLimit)/sizeof(KmerPosition<T>))+1);

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T>(par, subMat, seqDbr, totalKmersPerSplit, splits);
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T> *hashSeqPair = NULL;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = hashRanges.size();
    size_t fromSplit = 0;
    size_t splitCount = 1;
    mpiRank = MMseqsMPI::rank;
    // if split size is great than nodes than we have to
    // distribute all splits equally over all nodes
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(size_t i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }
    splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    for(size_t split = fromSplit; split < fromSplit+splitCount; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation<T>(totalKmers, hashRanges[split].first, hashRanges[split].second, splitFileName, seqDbr, par, subMat, true);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    for(size_t split = 0; split < hashRanges.size(); split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        Debug(Debug::INFO) << "Generate k-mers list for " << (split+1) <<" split\n";

        std::string splitFileNameDone = splitFileName + ".done";
        if(FileUtil::fileExists(splitFileNameDone.c_str()) == false){
            hashSeqPair = doComputation<T>(totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second, splitFileName, seqDbr, par, subMat, true);
        }

        splitFiles.push_back(splitFileName);
    }
#endif
    if(mpiRank == 0){
        std::vector<char> repSequence(seqDbr.getLastKey()+1);
        std::fill(repSequence.begin(), repSequence.end(), false);
        // write result
        DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed,
                     (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES );
        dbw.open();

        Timer timer;
        if(splits > 1) {
            seqDbr.unmapData();
            // Louis was here, added seqDbr as parameter to mergeKmerFilesAndOutput
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(dbw, splitFiles, repSequence);
            }else{
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(dbw, splitFiles, repSequence);
            }
            for(size_t i = 0; i < splitFiles.size(); i++){
                FileUtil::remove(splitFiles[i].c_str());
                std::string splitFilesDone = splitFiles[i] + ".done";
                FileUtil::remove(splitFilesDone.c_str());
            }
        } else {
            // Louis was here, added seqDbr as parameter to writeKmerMatcherResult
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1, true);
            }else{
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1, true);
            }
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";
        // add missing entries to the result (needed for clustering)

#pragma omp parallel num_threads(1)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                char buffer[100];
                unsigned int dbKey = seqDbr.getDbKey(id);
                if (repSequence[dbKey] == false) {
                    hit_t h;
                    h.prefScore = 0;
                    h.diagonal = 0;
                    h.seqId = dbKey;
                    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    // Louis was here
                    dbw.writeData(buffer, len, dbKey, thread_idx, seqDbr.getExtData(id));
                }
            }
        }
        dbw.close(false, false);
    }
    // free memory
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }

    return EXIT_SUCCESS;
}

template <int TYPE, typename T>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads, bool rymatch) {
    std::vector<size_t> threadOffsets;
    size_t splitSize = totalKmers/threads;
    threadOffsets.push_back(0);
    for(size_t thread = 1; thread < threads; thread++){
        size_t kmer = hashSeqPair[thread*splitSize].kmer;
        size_t repSeqId = static_cast<size_t>(kmer);
        repSeqId=BIT_SET(repSeqId, 63);
        bool wasSet = false;
        for(size_t pos = thread*splitSize; pos < totalKmers; pos++){
            size_t currSeqId = hashSeqPair[pos].kmer;
            currSeqId=BIT_SET(currSeqId, 63);
            if(repSeqId != currSeqId){
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }
        if(wasSet == false){
            threadOffsets.push_back(totalKmers - 1 );
        }
    }
    threadOffsets.push_back(totalKmers);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for(size_t thread = 0; thread < threads; thread++){
        std::string prefResultsOutString;
        prefResultsOutString.reserve(100000000);
        char buffer[100];
        size_t lastTargetId = SIZE_T_MAX;
        unsigned int writeSets = 0;
        size_t kmerPos=0;
        size_t repSeqId = SIZE_T_MAX;
        for(kmerPos = threadOffsets[thread]; kmerPos < threadOffsets[thread+1] && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++){
            size_t currKmer = hashSeqPair[kmerPos].kmer;
            int reverMask = 0;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                reverMask  = BIT_CHECK(currKmer, 63)==false;
                currKmer = BIT_CLEAR(currKmer, 63);
            }
            if(repSeqId != currKmer) {
                if (writeSets > 0) {
                    repSequence[repSeqId] = true;
                    dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), repSeqId, thread);
                }else{
                    if(repSeqId != SIZE_T_MAX) {
                        repSequence[repSeqId] = false;
                    }
                }
                lastTargetId = SIZE_T_MAX;
                prefResultsOutString.clear();
                repSeqId = currKmer;
                hit_t h;
                h.seqId = repSeqId;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                // TODO: error handling for len
                prefResultsOutString.append(buffer, len);
            }
            unsigned int targetId = hashSeqPair[kmerPos].id;
            T diagonal = hashSeqPair[kmerPos].pos;
            size_t kmerOffset = 0;
            T prevDiagonal = diagonal;
            size_t maxDiagonal = 0;
            size_t diagonalCnt = 0;
            size_t topScore =0;
            int bestReverMask = reverMask;
            // compute best diagonal and score for every group of target sequences
            while(lastTargetId != targetId
                  && kmerPos+kmerOffset < threadOffsets[thread+1]
                  && hashSeqPair[kmerPos+kmerOffset].id == targetId){
                if(prevDiagonal == hashSeqPair[kmerPos+kmerOffset].pos){
                    diagonalCnt++;
                }else{
                    diagonalCnt = 1;
                }
                if(diagonalCnt >= maxDiagonal){
                    diagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                    maxDiagonal = diagonalCnt;
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        bestReverMask = BIT_CHECK(hashSeqPair[kmerPos+kmerOffset].kmer, 63) == false;
                    }
                }
                prevDiagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                kmerOffset++;
                topScore++;
            }
            // remove similar double sequence hit
            if(targetId != repSeqId && lastTargetId != targetId ){
                ;
            }else{
                lastTargetId = targetId;
                continue;
            }
            hit_t h;
            h.seqId = targetId;
            h.prefScore = (bestReverMask) ? -topScore : topScore;
            h.diagonal = diagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
            lastTargetId = targetId;
            writeSets++;
        }
        if (writeSets > 0) {
            repSequence[repSeqId] = true;
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), repSeqId, thread);
        }else{
            if(repSeqId != SIZE_T_MAX) {
                repSequence[repSeqId] = false;
            }
        }
    }
}


int rymermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    int querySeqType = seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter *> *params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: " << seqDbr.getSize() << " type: " << seqDbr.getDbTypeName() << "\n";

    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        rymermatcherInner<short>(par, seqDbr);
    }
    else {
        rymermatcherInner<int>(par, seqDbr);
    }

    seqDbr.close();

    return EXIT_SUCCESS;
}

#undef SIZE_T_MAX
