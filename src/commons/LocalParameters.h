#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>
#include <MultiParam.h>
#include <algorithm>
#include <cfloat>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
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
    std::vector<MMseqsParameter *> ancientKmermatcherContigs;
    std::vector<MMseqsParameter *> ancientKmermatcherReads;


    int filterProteins;
    int deleteFilesInc;
    int minContigLen;
    int contigOutputMode;
    float clustSeqIdThr;
    float clustCovThr;
    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;
    bool dbMode;
    bool keepTarget;

    // ancient parameters
    float randomAlignPenal;
    float excessPenal;
    float likelihoodThreshold;
    float corrReadsRySeqId;
    float corrReadsSeqId;
    float corrContigSeqId;
    float mergeSeqIdThr;
    std::string ancientDamagePath;
    int numIterationsReads;
    bool ancientIncludeOnlyExtendReads;
    int ancientKmerSizeReads;
    bool ancientIncludeOnlyExtendContigs;
    int ancientKmerSizeContigs;
    float ancientKmersPerSequenceScale;
    int ancientkmersPerSequence;
    bool ancientUnsafe;
    float rySeqIdThr;

    MultiParam<int> multiNumIterations;
    MultiParam<int> multiKmerSize;
    MultiParam<int> multiAlnLenThr;
    MultiParam<float> multiSeqIdThr;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_MIN_CONTIG_LEN)
    PARAMETER(PARAM_CONTIG_OUTPUT_MODE)
    PARAMETER(PARAM_CLUST_MIN_SEQ_ID_THR)
    PARAMETER(PARAM_CLUST_C)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    PARAMETER(PARAM_MULTI_NUM_ITERATIONS)
    PARAMETER(PARAM_MULTI_K)
    PARAMETER(PARAM_MULTI_MIN_SEQ_ID)
    PARAMETER(PARAM_MULTI_MIN_ALN_LEN)
    PARAMETER(PARAM_DB_MODE)
    PARAMETER(PARAM_KEEP_TARGET)
    PARAMETER(PARAM_RAND_ALIGN)
    PARAMETER(PARAM_EXCESS_PENAL)
    PARAMETER(PARAM_CORR_THRESH)
    PARAMETER(PARAM_CORR_THRESH_SEQID)
    PARAMETER(PARAM_READ_EXT_THRESH)
    PARAMETER(PARAM_DAMAGE_PATH)
    PARAMETER(PARAM_NUM_ITERATIONS_READS)
    PARAMETER(PARAM_MIN_SEQ_MERGE_ID)
    PARAMETER(PARAM_CORR_CONTIG_MIN_SEQ_ID)
    PARAMETER(PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_READS)
    PARAMETER(PARAM_ANCIENT_PARAM_K_READS)
    PARAMETER(PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_CONTIGS)
    PARAMETER(PARAM_ANCIENT_PARAM_K_CONTIGS)
    PARAMETER(PARAM_ANCIENT_KMER_PER_SEQ_SCALE)
    PARAMETER(PARAM_ANCIENT_KMER_PER_SEQ)
    PARAMETER(PARAM_ANCIENT_MIN_RYSEQ_ID)
    PARAMETER(PARAM_ANCIENT_UNSAFE)


    // contig output
    static const int OUTPUT_ALL_CONTIGS = 0;
    static const int OUTPUT_ONLY_EXTENDED_CONTIGS = 1;
private:
    LocalParameters() :
            Parameters(),
            multiNumIterations(INT_MAX,INT_MAX),
            multiKmerSize(INT_MAX,INT_MAX),
            multiAlnLenThr(INT_MAX,INT_MAX),
            multiSeqIdThr(FLT_MAX,FLT_MAX),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "Delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
            PARAM_MIN_CONTIG_LEN(PARAM_MIN_CONTIG_LEN_ID, "--min-contig-len", "Minimum contig length", "Minimum length of assembled contig to output", typeid(int), (void *) &minContigLen, "^[1-9]{1}[0-9]*$"),
            PARAM_CONTIG_OUTPUT_MODE(PARAM_CONTIG_OUTPUT_MODE_ID, "--contig-output-mode", "Contig output mode", "Type of contigs:\n0: all\n1: only extended", typeid(int), (void *) &contigOutputMode, "^[0-1]{1}"),
            PARAM_CLUST_MIN_SEQ_ID_THR(PARAM_CLUST_MIN_SEQ_ID_THR_ID,"--clust-min-seq-id", "Clustering seq. id. threshold","Seq. id. threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)",typeid(float), (void *) &clustSeqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
            PARAM_CLUST_C(PARAM_CLUST_C_ID,"--clust-min-cov", "Clustering coverage threshold","Coverage threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)",typeid(float), (void *) &clustCovThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
            PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID,"--cycle-check", "Check for circular sequences", "Check for circular sequences (avoid over extension of circular or long repeated regions) ",typeid(bool), (void *) &cycleCheck, "", MMseqsParameter::COMMAND_MISC),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superfluous part of circular fragments (see --cycle-check)",typeid(bool), (void *) &chopCycle, "", MMseqsParameter::COMMAND_MISC),
            PARAM_MULTI_NUM_ITERATIONS(PARAM_MULTI_NUM_ITERATIONS_ID, "--num-iterations", "Number of assembly iterations","Number of assembly total iterations performed on nucleotide level (ignore protein level for ancient) (range 1-inf)",typeid(MultiParam<int>),(void *) &multiNumIterations, ""),
            PARAM_MULTI_K(PARAM_MULTI_K_ID, "-k", "k-mer length", "k-mer length (0: automatically set to optimum)", typeid(MultiParam<int>), (void *) &multiKmerSize, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
            PARAM_MULTI_MIN_SEQ_ID(PARAM_MULTI_MIN_SEQ_ID_ID, "--min-seq-id", "Seq. id. threshold", "Overlap sequence identity threshold [0.0, 1.0]", typeid(MultiParam<float>), (void *) &multiSeqIdThr, "", MMseqsParameter::COMMAND_ALIGN),
            PARAM_MULTI_MIN_ALN_LEN(PARAM_MULTI_MIN_ALN_LEN_ID, "--min-aln-len", "Min alignment length", "Minimum alignment length (range 0-INT_MAX)", typeid(MultiParam<int>), (void *) &multiAlnLenThr, "", MMseqsParameter::COMMAND_ALIGN),
            PARAM_DB_MODE(PARAM_DB_MODE_ID, "--db-mode", "Input is database", "Input is database", typeid(bool), (void *) &dbMode, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_KEEP_TARGET(PARAM_KEEP_TARGET_ID, "--keep-target", "Keep target sequences for the next iteration", "Keep target sequences", typeid(bool), (void*) &keepTarget, "", MMseqsParameter::COMMAND_MISC), 
            PARAM_RAND_ALIGN(PARAM_RAND_ALIGN_ID, "--ext-random-align", "Random alignment (ancient)", "Use either: 0.8 or 0.9 (ancient)", typeid(float), (void *) &randomAlignPenal, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_EXCESS_PENAL(PARAM_EXCESS_PENAL_ID, "--excess-penalty", "Penalize short overlaps (ancient)", "Use float: 0.25 to 0.5 (ancient)", typeid(float), (void *) &excessPenal, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_CORR_THRESH(PARAM_CORR_THRESH_ID, "--min-ryseq-id-corr-reads", "Seq. ident. in RY-space to increase precision of correction (ancient)" , "Min. RY-mer space seq. ident in correction phase. Range 0-1 (ancient)", typeid(float), (void *) &corrReadsRySeqId, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_CORR_THRESH_SEQID(PARAM_CORR_THRESH_SEQID_ID, "--min-seqid-corr-reads", "Seq. ident. to increase precision of correction (ancient)" , "Min. seq. ident. in correction phase. Range 0-1 (ancient)", typeid(float), (void *) &corrReadsSeqId, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_READ_EXT_THRESH(PARAM_READ_EXT_THRESH_ID, "--likelihood-ratio-threshold", "Min. odds ratio threshold for read extension (ancient)" , " Min. odds ratio to accept read extension. Range 0-1 (ancient)", typeid(float), (void *) &likelihoodThreshold, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_DAMAGE_PATH(PARAM_DAMAGE_PATH_ID, "--ancient-damage", "Path to deamination patterns of ancient DNA (ancient)", "Path to damage matrix (ancient)", typeid(std::string), (void *) &ancientDamagePath, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_NUM_ITERATIONS_READS(PARAM_NUM_ITERATIONS_READS_ID, "--num-iter-reads-only", "Number of assembly iterations with raw reads only (ancient)", "Raw reads only: Number of assembly iterations performed on nucleotide level (ancient)", typeid(int), (void *)&numIterationsReads, "^[1-9]{1}[0-9]*$"),
            PARAM_MIN_SEQ_MERGE_ID(PARAM_MIN_SEQ_MERGE_ID_ID, "--min-merge-seq-id", "Seq. id. threshold during merging (ancient)", "Min. seq. ident. of contig overlaps (ancient) (range 0.0-1.0)", typeid(float), (void *) &mergeSeqIdThr, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_CORR_CONTIG_MIN_SEQ_ID(PARAM_CORR_CONTIG_MIN_SEQ_ID_ID, "--min-seqid-corr-contigs", "Seq. id. threshold during contig correction (ancient)", "Min. seq. ident. for contig correction (ancient) (range 0.0-1.0)", typeid(float), (void *) &corrContigSeqId, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_READS(PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_READS_ID, "--include-only-extendable-ancient-reads", "Include only extendable in read extension (ancient)", "Include only extendable (reads onl, ancient)", typeid(bool), (void *) &ancientIncludeOnlyExtendReads, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_PARAM_K_READS(PARAM_ANCIENT_PARAM_K_READS_ID, "--k-ancient-reads", "k-mer length reads (ancient)", "k-mer length read step (ancient)", typeid(int), (void *) &ancientKmerSizeReads, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_CONTIGS(PARAM_ANCIENT_INCLUDE_ONLY_EXTENDABLE_CONTIGS_ID, "--include-only-extendable-ancient-contigs", "Include only extendable in contig merging (ancient)", "Include only extendable (contigs, ancient)", typeid(bool), (void *) &ancientIncludeOnlyExtendContigs, "", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_PARAM_K_CONTIGS(PARAM_ANCIENT_PARAM_K_CONTIGS_ID, "--k-ancient-contigs", "k-mer length contigs (ancient)", "k-mer length contig step (ancient)", typeid(int), (void *) &ancientKmerSizeContigs, "^[0-9]{1}[0-9]*$", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_KMER_PER_SEQ_SCALE(PARAM_ANCIENT_KMER_PER_SEQ_SCALE_ID, "--kmer-per-seq-scale-ancient", "Scale k-mers per sequence (ancient)", "Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen (ancient)", typeid(MultiParam<float>), (void *) &kmersPerSequenceScale, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_KMER_PER_SEQ(PARAM_ANCIENT_KMER_PER_SEQ_ID, "--kmer-per-seq-ancient", "k-mers per sequence (ancient)", "k-mers per sequence (ancient)", typeid(int), (void *) &kmersPerSequence, "^[1-9]{1}[0-9]*$", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_MIN_RYSEQ_ID(PARAM_ANCIENT_MIN_RYSEQ_ID_ID, "--min-ryseq-id", "RY-Seq. id. threshold (ancient)", "List matches above this sequence identity in RY-mer space (ancient) (range 0.0-1.0)", typeid(float), (void *) &rySeqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_EXPERT),
            PARAM_ANCIENT_UNSAFE(PARAM_ANCIENT_UNSAFE_ID, "--unsafe", "Maximize contig length (safe vs. unsafe mode)", "Maximize the contig length, but higher misassembly rate (ancient) ", typeid(bool), (void *) &ancientUnsafe, "", MMseqsParameter::COMMAND_EXPERT)
            {

        // assembleresult
        assembleresults.push_back(&PARAM_MIN_SEQ_ID);
        assembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        assembleresults.push_back(&PARAM_KEEP_TARGET);
        assembleresults.push_back(&PARAM_THREADS);
        assembleresults.push_back(&PARAM_V);
        assembleresults.push_back(
                &PARAM_RESCORE_MODE); //temporary added until assemble and nuclassemble use same rescoremode
        assembleresults.push_back(&PARAM_RAND_ALIGN);
        assembleresults.push_back(&PARAM_EXCESS_PENAL);
        assembleresults.push_back(&PARAM_CORR_THRESH);
        assembleresults.push_back(&PARAM_CORR_THRESH_SEQID);
        assembleresults.push_back(&PARAM_READ_EXT_THRESH);
        assembleresults.push_back(&PARAM_MIN_SEQ_MERGE_ID);
        assembleresults.push_back(&PARAM_CORR_CONTIG_MIN_SEQ_ID);
        assembleresults.push_back(&PARAM_DAMAGE_PATH);
        assembleresults.push_back(&PARAM_ANCIENT_UNSAFE);


        extractorfssubset.push_back(&PARAM_TRANSLATION_TABLE);
        extractorfssubset.push_back(&PARAM_USE_ALL_TABLE_STARTS);
        extractorfssubset.push_back(&PARAM_THREADS);
        extractorfssubset.push_back(&PARAM_V);

        filternoncoding.push_back(&PARAM_PROTEIN_FILTER_THRESHOLD);
        filternoncoding.push_back(&PARAM_THREADS);
        filternoncoding.push_back(&PARAM_V);

        //cyclecheck
        cyclecheck.push_back(&PARAM_MAX_SEQ_LEN);
        cyclecheck.push_back(&PARAM_CHOP_CYCLE);
        cyclecheck.push_back(&PARAM_THREADS);
        cyclecheck.push_back(&PARAM_V);

        //createhdb
        createhdb.push_back(&PARAM_COMPRESSED);
        createhdb.push_back(&PARAM_V);

        //reduceredundancy (subset of clustering parameters which have to be adjusted)
        reduceredundancy.push_back(&PARAM_ALPH_SIZE);
        reduceredundancy.push_back(&PARAM_CLUSTER_MODE);
        reduceredundancy.push_back(&PARAM_K);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ_SCALE);
        reduceredundancy.push_back(&PARAM_IGNORE_MULTI_KMER);
        reduceredundancy.push_back(&PARAM_MIN_SEQ_ID);
        reduceredundancy.push_back(&PARAM_COV_MODE);
        reduceredundancy.push_back(&PARAM_C);
        reduceredundancy.push_back(&PARAM_MAX_SEQ_LEN);
        reduceredundancy.push_back(&PARAM_WRAPPED_SCORING);
        reduceredundancy.push_back(&PARAM_GAP_OPEN);
        reduceredundancy.push_back(&PARAM_GAP_EXTEND);
        reduceredundancy.push_back(&PARAM_ZDROP);
        reduceredundancy.push_back(&PARAM_THREADS);
        reduceredundancy.push_back(&PARAM_REMOVE_TMP_FILES);

        // kmermatcher contigs (subset of kmermatcher parameters which have to be adjusted)
        ancientKmermatcherContigs.push_back(&PARAM_INCLUDE_ONLY_EXTENDABLE);
        ancientKmermatcherContigs.push_back(&PARAM_K);
        ancientKmermatcherContigs.push_back(&PARAM_KMER_PER_SEQ);
        ancientKmermatcherContigs.push_back(&PARAM_KMER_PER_SEQ_SCALE);
        ancientKmermatcherContigs = combineList(kmermatcher,ancientKmermatcherContigs);

        // kmermatcher contigs (subset of kmermatcher parameters which have to be adjusted)
        ancientKmermatcherReads.push_back(&PARAM_INCLUDE_ONLY_EXTENDABLE);
        ancientKmermatcherReads.push_back(&PARAM_K);
        ancientKmermatcherReads.push_back(&PARAM_KMER_PER_SEQ);
        ancientKmermatcherReads.push_back(&PARAM_KMER_PER_SEQ_SCALE);
        ancientKmermatcherReads = combineList(kmermatcher,ancientKmermatcherReads);

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
        nuclassembleworkflow = combineList(nuclassembleworkflow, rescorediagonal);
        nuclassembleworkflow = combineList(nuclassembleworkflow, assembleresults);
        nuclassembleworkflow = combineList(nuclassembleworkflow, cyclecheck);
        nuclassembleworkflow = combineList(nuclassembleworkflow, ancientKmermatcherContigs);

        nuclassembleworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassembleworkflow.push_back(&PARAM_MIN_CONTIG_LEN);
        nuclassembleworkflow.push_back(&PARAM_CONTIG_OUTPUT_MODE);
        nuclassembleworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassembleworkflow.push_back(&PARAM_DB_MODE);
        nuclassembleworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassembleworkflow.push_back(&PARAM_DELETE_TMP_INC);
        nuclassembleworkflow.push_back(&PARAM_RUNNER);
        nuclassembleworkflow.push_back(&PARAM_NUM_ITERATIONS_READS);
        nuclassembleworkflow.push_back(&PARAM_ANCIENT_MIN_RYSEQ_ID);
        nuclassembleworkflow.push_back(&PARAM_ANCIENT_PARAM_K_READS);
        nuclassembleworkflow.push_back(&PARAM_ANCIENT_PARAM_K_CONTIGS);

        // guidedassembleresults
        guidedassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        guidedassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        guidedassembleresults.push_back(&PARAM_KEEP_TARGET);
        guidedassembleresults.push_back(&PARAM_RESCORE_MODE);
        guidedassembleresults.push_back(&PARAM_THREADS);
        guidedassembleresults.push_back(&PARAM_V);

        // guidedNuclAssembleworkflow
        guidedNuclAssembleworkflow = combineList(extractorfs, guidedassembleresults);
        guidedNuclAssembleworkflow = combineList(guidedNuclAssembleworkflow, nuclassembleworkflow);
        guidedNuclAssembleworkflow = combineList(guidedNuclAssembleworkflow, reduceredundancy);
        guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_MIN_SEQ_ID_THR);
        guidedNuclAssembleworkflow.push_back(&PARAM_CLUST_C);

        // guidedNuclAssembleworkflow special parameter: replace with MultiParam to make aa and nucl values independent
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_K);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_K);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_MIN_SEQ_ID);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_MIN_SEQ_ID);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_MIN_ALN_LEN);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_MIN_ALN_LEN);
        guidedNuclAssembleworkflow = removeParameter(guidedNuclAssembleworkflow, PARAM_NUM_ITERATIONS);
        guidedNuclAssembleworkflow.push_back(&PARAM_MULTI_NUM_ITERATIONS);

        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;
        clustSeqIdThr = 0.97;
        clustCovThr = 0.99;
        minContigLen = 500;
        contigOutputMode = 1;
        chopCycle = true;
        cycleCheck = true;
        dbMode = false;
        keepTarget = true;

        randomAlignPenal = 0.85;
        excessPenal = 0.0625;
        likelihoodThreshold = 0.5;
        rySeqIdThr = 0.99;
        corrReadsRySeqId = 0.99;
        corrReadsSeqId = 0.9;
        corrContigSeqId = 0.9;
        mergeSeqIdThr = 0.99;
        ancientDamagePath = "";
        numIterationsReads = 5;

        ancientIncludeOnlyExtendReads = false;
        ancientIncludeOnlyExtendContigs = true;
        ancientKmerSizeReads = 20;
        ancientKmerSizeContigs = 22;
        ancientKmersPerSequenceScale = 0.2;
        ancientkmersPerSequence = 200;
        ancientUnsafe = false;

        multiNumIterations = MultiParam<int>(12, 1);
        multiKmerSize = MultiParam<int>(14, 20);
        multiAlnLenThr = MultiParam<int>(0, 0);
        multiSeqIdThr = MultiParam<float>(0.97, 0.97);

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
