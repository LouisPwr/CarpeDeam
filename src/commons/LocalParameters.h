#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <MultiParam.h>
#include <Parameters.h>

#include <algorithm>
#include <cfloat>

class LocalParameters : public Parameters {
   public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters &getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters &>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter *> assembleworkflow;
    std::vector<MMseqsParameter *> nuclassembleworkflow;
    std::vector<MMseqsParameter *> guidedNuclAssembleworkflow;

    std::vector<MMseqsParameter *> assembleresults;
    std::vector<MMseqsParameter *> cyclecheck;
    std::vector<MMseqsParameter *> createhdb;
    std::vector<MMseqsParameter *> extractorfssubset;
    std::vector<MMseqsParameter *> filternoncoding;
    std::vector<MMseqsParameter *> guidedassembleresults;
    std::vector<MMseqsParameter *> reduceredundancy;

    int filterProteins;
    int deleteFilesInc;
    int minContigLen;
    int usePrefilter;
    int prefilterNumIterations;
    float clustSeqIdThr;
    float clustRySeqIdThr;
    // int clustRySize;
    float clustCovThr;

    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;
    bool dbMode;
    int prefilterSpacedKmer;
    std::string prefilterSpacedKmerPattern;
    int prefilterExactKmerMatching;
    int prefilterMaskMode;
    int prefilterCompBiasCorrection;
    float prefilterSensitivity;
    int prefilterKmerSize;
    size_t prefilterMaxResListLen;                // Maximal result list length per query

    float randomAlignPenal;
    float correctionThreshold;
    float mergeSeqIdThr;
    std::string ancientDamagePath;


    MultiParam<char*> prefilterScoringMatrixFile;       // path to scoring matrix
    MultiParam<int> multiNumIterations;
    MultiParam<int> multiKmerSize;
    MultiParam<int> multiAlnLenThr;
    MultiParam<float> multiSeqIdThr;
    MultiParam<int> multiSpacedKmer;
    MultiParam<char *> multiSpacedKmerPattern;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_MIN_CONTIG_LEN)
    PARAMETER(PARAM_CLUST_MIN_SEQ_ID_THR)
    PARAMETER(PARAM_CLUST_MIN_RYSEQ_ID_THR)
    // PARAMETER(PARAM_CLUST_RYMER_LEN)
    PARAMETER(PARAM_CLUST_C)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    PARAMETER(PARAM_MULTI_NUM_ITERATIONS)
    PARAMETER(PARAM_MULTI_K)
    PARAMETER(PARAM_MULTI_MIN_SEQ_ID)
    PARAMETER(PARAM_MULTI_MIN_ALN_LEN)
    PARAMETER(PARAM_DB_MODE)
    PARAMETER(PARAM_USE_PREFILTER)
    PARAMETER(PARAM_PREFILTER_NUM_ITERATIONS)
    PARAMETER(PARAM_PREFILTER_SUB_MAT)
    PARAMETER(PARAM_PREFILTER_S)
    PARAMETER(PARAM_PREFILTER_K)
    PARAMETER(PARAM_PREFILTER_SPACED_KMER_MODE)
    PARAMETER(PARAM_PREFILTER_SPACED_KMER_PATTERN)
    PARAMETER(PARAM_PREFILTER_EXACT_KMER_MATCHING)
    PARAMETER(PARAM_PREFILTER_MASK_RESIDUES)
    PARAMETER(PARAM_PREFILTER_NO_COMP_BIAS_CORR)
    PARAMETER(PARAM_PREFILTER_MAX_SEQS)
    PARAMETER(PARAM_MULTI_SPACED_KMER_MODE)
    PARAMETER(PARAM_MULTI_SPACED_KMER_PATTERN)
    PARAMETER(PARAM_RAND_ALIGN)
    PARAMETER(PARAM_CORR_THRESH)
    PARAMETER(PARAM_DAMAGE_PATH)
    PARAMETER(PARAM_MIN_SEQ_MERGE_ID)

   private:
    LocalParameters() : Parameters(),
                        prefilterScoringMatrixFile("INVALID", "INVALID"),
                        multiNumIterations(INT_MAX, INT_MAX),
                        multiKmerSize(INT_MAX, INT_MAX),
                        multiAlnLenThr(INT_MAX, INT_MAX),
                        multiSeqIdThr(FLT_MAX, FLT_MAX),
                        multiSpacedKmer(INT_MAX, INT_MAX),
                        multiSpacedKmerPattern("", ""),
                        PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID, "--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]", typeid(int), (void *)&filterProteins, "^[0-1]{1}$"),
                        PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID, "--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]", typeid(float), (void *)&proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
                        PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID, "--delete-tmp-inc", "Delete temporary files incremental", "Delete temporary files incremental [0,1]", typeid(int), (void *)&deleteFilesInc, "^[0-1]{1}$", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
                        PARAM_MIN_CONTIG_LEN(PARAM_MIN_CONTIG_LEN_ID, "--min-contig-len", "Minimum contig length", "Minimum length of assembled contig to output", typeid(int), (void *)&minContigLen, "^[1-9]{1}[0-9]*$"),
                        PARAM_CLUST_MIN_SEQ_ID_THR(PARAM_CLUST_MIN_SEQ_ID_THR_ID, "--clust-min-seq-id", "Clustering seq. id. threshold", "Seq. id. threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)", typeid(float), (void *)&clustSeqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
                        PARAM_CLUST_MIN_RYSEQ_ID_THR(PARAM_CLUST_MIN_RYSEQ_ID_THR_ID, "--clust-min-ry-seq-id", "Clustering seq. id. threshold in RY-space", "RY-Seq. id. threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)", typeid(float), (void *)&clustRySeqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),      
                        // PARAM_CLUST_RYMER_LEN(PARAM_CLUST_RYMER_LEN_ID, "--clust-rymer-len", "Lenght of RY-mers used for clustering", "RY-mer length for redundancy reduction", typeid(int), (void *)&clustRySize, "", MMseqsParameter::COMMAND_CLUST),
                        PARAM_CLUST_C(PARAM_CLUST_C_ID, "--clust-min-cov", "Clustering coverage threshold", "Coverage threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)", typeid(float), (void *)&clustCovThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
                        PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID, "--cycle-check", "Check for circular sequences", "Check for circular sequences (avoid over extension of circular or long repeated regions) ", typeid(bool), (void *)&cycleCheck, "", MMseqsParameter::COMMAND_MISC),
                        PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID, "--chop-cycle", "Chop Cycle", "Remove superfluous part of circular fragments (see --cycle-check)", typeid(bool), (void *)&chopCycle, "", MMseqsParameter::COMMAND_MISC),
                        PARAM_MULTI_NUM_ITERATIONS(PARAM_MULTI_NUM_ITERATIONS_ID, "--num-iterations", "Number of assembly iterations", "Number of assembly iterations performed on nucleotide level,protein level (range 1-inf)", typeid(MultiParam<int>), (void *)&multiNumIterations, ""),
                        PARAM_MULTI_K(PARAM_MULTI_K_ID, "-k", "k-mer length", "k-mer length (0: automatically set to optimum)", typeid(MultiParam<int>), (void *)&multiKmerSize, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
                        PARAM_MULTI_MIN_SEQ_ID(PARAM_MULTI_MIN_SEQ_ID_ID, "--min-seq-id", "Seq. id. threshold", "Overlap sequence identity threshold [0.0, 1.0]", typeid(MultiParam<float>), (void *)&multiSeqIdThr, "", MMseqsParameter::COMMAND_ALIGN),
                        PARAM_MULTI_MIN_ALN_LEN(PARAM_MULTI_MIN_ALN_LEN_ID, "--min-aln-len", "Min alignment length", "Minimum alignment length (range 0-INT_MAX)", typeid(MultiParam<int>), (void *)&multiAlnLenThr, "", MMseqsParameter::COMMAND_ALIGN),
                        PARAM_DB_MODE(PARAM_DB_MODE_ID, "--db-mode", "Input is database", "Input is database", typeid(bool), (void *)&dbMode, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_USE_PREFILTER(PARAM_USE_PREFILTER_ID, "--use-prefilter", "Use prefilter for assembly", "Use prefilter to assemble super short reads [0,1]", typeid(int), (void *)&usePrefilter, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_NUM_ITERATIONS(PARAM_PREFILTER_NUM_ITERATIONS_ID, "--prefilter-iterations", "Number of assemly iteration using the prefilter", "Number of assembly iterations using the prefilter before changing to kmermatcher", typeid(int), (void *)&prefilterNumIterations, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_SUB_MAT(PARAM_PREFILTER_SUB_MAT_ID, "--prefilter-sub-mat", "Prefilter substitution matrix", "Substitution matrix file", typeid(MultiParam<char*>), (void *) &prefilterScoringMatrixFile, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_S(PARAM_PREFILTER_S_ID, "--prefilter-s", "Prefilter sensitivity", "Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive", typeid(float), (void *) &prefilterSensitivity, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_K(PARAM_PREFILTER_K_ID, "--prefilter-k", "Prefilter k-mer length", "Prefilter k-mer length (0: automatically set to optimum)", typeid(int), (void *)&prefilterKmerSize, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_SPACED_KMER_MODE(PARAM_PREFILTER_SPACED_KMER_MODE_ID, "--prefilter-spaced-kmer-mode", "Prefilter spaced k-mers", "0: use consecutive positions in k-mers; 1: use spaced k-mers", typeid(int), (void *) &prefilterSpacedKmer, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_SPACED_KMER_PATTERN(PARAM_PREFILTER_SPACED_KMER_PATTERN_ID, "--prefilter-spaced-kmer-pattern", "Prefilter user-specified spaced k-mer pattern", "User-specified spaced k-mer pattern", typeid(std::string), (void *) &prefilterSpacedKmerPattern, "^1[01]*1$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_EXACT_KMER_MATCHING(PARAM_PREFILTER_EXACT_KMER_MATCHING_ID, "--prefilter-exact-kmer-matching", "Prefilter exact k-mer matching", "Extract only exact k-mers for matching (range 0-1)", typeid(int), (void *) &prefilterExactKmerMatching, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_MASK_RESIDUES(PARAM_PREFILTER_MASK_RESIDUES_ID, "--prefilter-mask", "Prefilter mask residues", "Mask sequences in k-mer stage: 0: w/o low complexity masking, 1: with low complexity masking", typeid(int), (void *) &prefilterMaskMode, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_NO_COMP_BIAS_CORR(PARAM_PREFILTER_NO_COMP_BIAS_CORR_ID, "--prefilter-comp-bias-corr", "Prefilter compositional bias", "Correct for locally biased amino acid composition (range 0-1)", typeid(int), (void *) &prefilterCompBiasCorrection, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_PREFILTER_MAX_SEQS(PARAM_PREFILTER_MAX_SEQS_ID, "--prefilter-max-seqs", "Prefilter max results per query", "Maximum results per query sequence allowed to pass the prefilter (affects sensitivity)", typeid(size_t), (void *) &prefilterMaxResListLen, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_MULTI_SPACED_KMER_MODE(PARAM_MULTI_SPACED_KMER_MODE_ID, "--spaced-kmer-mode", "Spaced k-mers", "0: use consecutive positions in k-mers; 1: use spaced k-mers", typeid(MultiParam<int>), (void *)&multiSpacedKmer, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_MULTI_SPACED_KMER_PATTERN(PARAM_MULTI_SPACED_KMER_PATTERN_ID, "--spaced-kmer-pattern", "User-specified spaced k-mer pattern", "User-specified spaced k-mer pattern", typeid(MultiParam<char *>), (void *)&multiSpacedKmerPattern, "", MMseqsParameter::COMMAND_EXPERT),

                        PARAM_RAND_ALIGN(PARAM_RAND_ALIGN_ID, "--ext-random-align", "random alignment", "Use either: 0.5, 0.75, 0.625", typeid(float), (void *) &randomAlignPenal, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_CORR_THRESH(PARAM_CORR_THRESH_ID, "--correction-min-ry-seqid", "Seq. ident. in RY-space to increase precision of correction" , "Range 0-1", typeid(float), (void *) &correctionThreshold, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_DAMAGE_PATH(PARAM_DAMAGE_PATH_ID, "--ancient-damage", "path to deamination patterns of ancient DNA", "Path to dir", typeid(std::string), (void *) &ancientDamagePath, "", MMseqsParameter::COMMAND_EXPERT),
                        PARAM_MIN_SEQ_MERGE_ID(PARAM_MIN_SEQ_MERGE_ID_ID, "--min-merge-seq-id", "Seq. id. threshold during merging", "List matches above this sequence identity (for merging contigs) (range 0.0-1.0)", typeid(float), (void *) &mergeSeqIdThr, "", MMseqsParameter::COMMAND_EXPERT){

        // assembleresult
        assembleresults.push_back(&PARAM_MIN_SEQ_ID);
        assembleresults.push_back(&PARAM_MIN_RYSEQ_ID);
        assembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        assembleresults.push_back(&PARAM_THREADS);
        assembleresults.push_back(&PARAM_V);
        assembleresults.push_back(
            &PARAM_RESCORE_MODE);  // temporary added until assemble and nuclassemble use same rescoremode
        assembleresults.push_back(&PARAM_RAND_ALIGN);
        assembleresults.push_back(&PARAM_CORR_THRESH);
        assembleresults.push_back(&PARAM_MIN_SEQ_MERGE_ID);
        assembleresults.push_back(&PARAM_DAMAGE_PATH);

        extractorfssubset.push_back(&PARAM_TRANSLATION_TABLE);
        extractorfssubset.push_back(&PARAM_USE_ALL_TABLE_STARTS);
        extractorfssubset.push_back(&PARAM_THREADS);
        extractorfssubset.push_back(&PARAM_V);

        filternoncoding.push_back(&PARAM_PROTEIN_FILTER_THRESHOLD);
        filternoncoding.push_back(&PARAM_THREADS);
        filternoncoding.push_back(&PARAM_V);

        // cyclecheck
        cyclecheck.push_back(&PARAM_MAX_SEQ_LEN);
        cyclecheck.push_back(&PARAM_CHOP_CYCLE);
        cyclecheck.push_back(&PARAM_THREADS);
        cyclecheck.push_back(&PARAM_V);

        // createhdb
        createhdb.push_back(&PARAM_COMPRESSED);
        createhdb.push_back(&PARAM_V);

        // reduceredundancy (subset of clustering parameters which have to be adjusted)
        reduceredundancy.push_back(&PARAM_ALPH_SIZE);
        reduceredundancy.push_back(&PARAM_CLUSTER_MODE);
        reduceredundancy.push_back(&PARAM_K);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ_SCALE);
        reduceredundancy.push_back(&PARAM_IGNORE_MULTI_KMER);
        reduceredundancy.push_back(&PARAM_MIN_SEQ_ID);
        reduceredundancy.push_back(&PARAM_MIN_RYSEQ_ID);
        reduceredundancy.push_back(&PARAM_COV_MODE);
        reduceredundancy.push_back(&PARAM_C);
        reduceredundancy.push_back(&PARAM_MAX_SEQ_LEN);
        reduceredundancy.push_back(&PARAM_WRAPPED_SCORING);
        reduceredundancy.push_back(&PARAM_GAP_OPEN);
        reduceredundancy.push_back(&PARAM_GAP_EXTEND);
        reduceredundancy.push_back(&PARAM_ZDROP);
        reduceredundancy.push_back(&PARAM_THREADS);
        reduceredundancy.push_back(&PARAM_REMOVE_TMP_FILES);

        // assembler workflow
        assembleworkflow = combineList(createdb, kmermatcher);
        assembleworkflow = combineList(assembleworkflow, rescorediagonal);
        assembleworkflow = combineList(assembleworkflow, extractorfs);
        assembleworkflow = combineList(assembleworkflow, assembleresults);
        assembleworkflow = combineList(assembleworkflow, filternoncoding);

        assembleworkflow.push_back(&PARAM_FILTER_PROTEINS);
        assembleworkflow.push_back(&PARAM_NUM_ITERATIONS);
        assembleworkflow.push_back(&PARAM_DELETE_TMP_INC);
        assembleworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        assembleworkflow.push_back(&PARAM_RUNNER);

        // nuclassembler workflow
        nuclassembleworkflow = combineList(createdb, kmermatcher);
        nuclassembleworkflow = combineList(nuclassembleworkflow, prefilter);
        nuclassembleworkflow = combineList(nuclassembleworkflow, rescorediagonal);
        nuclassembleworkflow = combineList(nuclassembleworkflow, assembleresults);
        nuclassembleworkflow = combineList(nuclassembleworkflow, cyclecheck);

        nuclassembleworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassembleworkflow.push_back(&PARAM_MIN_CONTIG_LEN);
        nuclassembleworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassembleworkflow.push_back(&PARAM_DB_MODE);
        nuclassembleworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassembleworkflow.push_back(&PARAM_DELETE_TMP_INC);
        nuclassembleworkflow.push_back(&PARAM_RUNNER);

        // guidedassembleresults
        guidedassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        guidedassembleresults.push_back(&PARAM_MIN_RYSEQ_ID);
        guidedassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        guidedassembleresults.push_back(&PARAM_RESCORE_MODE);
        guidedassembleresults.push_back(&PARAM_THREADS);
        guidedassembleresults.push_back(&PARAM_V);

        // guidedNuclAssembleworkflow
        guidedNuclAssembleworkflow = combineList(extractorfs, guidedassembleresults);
        guidedNuclAssembleworkflow = combineList(guidedNuclAssembleworkflow, nuclassembleworkflow);
        guidedNuclAssembleworkflow = combineList(guidedNuclAssembleworkflow, reduceredundancy);
        guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_MIN_SEQ_ID_THR);
        guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_MIN_RYSEQ_ID_THR);
        // guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_RYMER_LEN);
        guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_C);
        guidedNuclAssembleworkflow.push_back(&PARAM_USE_PREFILTER);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_NUM_ITERATIONS);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SUB_MAT);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_S);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_K);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SPACED_KMER_MODE);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SPACED_KMER_PATTERN);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_EXACT_KMER_MATCHING);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_MASK_RESIDUES);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_NO_COMP_BIAS_CORR);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_MAX_SEQS);

        // guidedNuclAssembleworkflow special parameter: replace with MultiParam to make aa and nucl values independent
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_K);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_K);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_MIN_SEQ_ID);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_MIN_SEQ_ID);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_MIN_ALN_LEN);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_MIN_ALN_LEN);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_NUM_ITERATIONS);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_NUM_ITERATIONS);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_SUB_MAT);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SUB_MAT);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_S);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_S);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_K);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_K);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_SPACED_KMER_MODE);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SPACED_KMER_MODE);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_SPACED_KMER_PATTERN);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_SPACED_KMER_PATTERN);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_SPACED_KMER_MODE);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_SPACED_KMER_MODE);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_SPACED_KMER_PATTERN);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_SPACED_KMER_PATTERN);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_EXACT_KMER_MATCHING);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_EXACT_KMER_MATCHING);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_MASK_RESIDUES);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_MASK_RESIDUES);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_NO_COMP_BIAS_CORR);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_NO_COMP_BIAS_CORR);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_PREFILTER_MAX_SEQS);
        guidedNuclAssembleworkflow.push_back(&PARAM_PREFILTER_MAX_SEQS);


        filterProteins = 1;
        deleteFilesInc = 1;
        usePrefilter = 0;
        proteinFilterThreshold = 0.2;
        clustSeqIdThr = 0.97;
        clustRySeqIdThr = 0.99;
        // clustRySize = 32;
        clustCovThr = 0.99;
        minContigLen = 1000;
        chopCycle = true;
        cycleCheck = true;
        dbMode = false;
        prefilterNumIterations = 3;
        prefilterSpacedKmer = 0;
        prefilterSensitivity = 4.0;
        prefilterKmerSize = 14;
        prefilterMaxResListLen = 300;
        prefilterExactKmerMatching = false;
        prefilterMaskMode = false;
        prefilterCompBiasCorrection = false;

        multiNumIterations = MultiParam<int>(8, 4);
        multiKmerSize = MultiParam<int>(14, 22);
        multiAlnLenThr = MultiParam<int>(0, 0);
        multiSeqIdThr = MultiParam<float>(0.97, 0.97);
        multiSpacedKmer = MultiParam<int>(0, 0);
        prefilterScoringMatrixFile =  MultiParam<char*>("blosum62.out", "nucleotide.out");

        randomAlignPenal = 0.675;
        correctionThreshold = 0.99;
        mergeSeqIdThr = 0.99;
        ancientDamagePath = "";

    }
    LocalParameters(LocalParameters const &);
    ~LocalParameters(){};
    void operator=(LocalParameters const &);
};

#endif
