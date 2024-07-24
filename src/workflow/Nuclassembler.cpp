#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "nuclassemble.sh.h"

void setNuclAssemblerWorkflowDefaults(LocalParameters *p) {

    p->numIterations = 10;
    p->numIterationsReads = 4;
    p->kmerSize = 22;
    p->seqIdThr = 0.9;
    p->mergeSeqIdThr = 0.99;
    p->alphabetSize = 5;

    p->covThr = 0.0;
    p->evalThr = 0.001;
    p->maskMode = 0;
    p->kmersPerSequence = 200;
    p->kmersPerSequenceScale = 0.2;
    p->spacedKmer = false;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = false;
    p->addBacktrace = false;
    p->rescoreMode = Parameters::RESCORE_MODE_END_TO_END_ALIGNMENT;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->maxSeqLen = 300000;
    p->cycleCheck = true;
    p->chopCycle = true;

}

int nuclassemble(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setNuclAssemblerWorkflowDefaults(&par);

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold (range 0.0-1.0)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations (range 1-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below (range 0.0-inf)", NULL, 0);

    // make parameter visible in short help
    par.overrideParameterDescription( par.PARAM_MAX_SEQ_LEN, NULL, NULL, MMseqsParameter::COMMAND_COMMON);


    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SORT_RESULTS.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    CommandCaller cmd;
    if(par.filenames.size() < 3) {
        Debug(Debug::ERROR) << "Too few input files provided.\n";
        return EXIT_FAILURE;
    }
    else if ((par.filenames.size() - 2) % 2 == 0) {
        cmd.addVariable("PAIRED_END", "1"); // paired end reads
    } else {
        if (par.filenames.size() != 3) {
            Debug(Debug::ERROR) << "Too many input files provided.\n";
            Debug(Debug::ERROR) << "For paired-end input provide READSETA_1.fastq READSETA_2.fastq ... OUTPUT.fasta tmpDir\n";
            Debug(Debug::ERROR) << "For single input use READSET.fast(q|a) OUTPUT.fasta tmpDir\n";
            return EXIT_FAILURE;
        }
        cmd.addVariable("PAIRED_END", NULL); // single end reads
    }

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);

    char *p = realpath(tmpDir.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get real path of " << tmpDir << "!\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("TMP_PATH", p);
    par.filenames.pop_back();
    free(p);

    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_INCREMENTAL_TMP", par.deleteFilesInc ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    cmd.addVariable("NUM_IT_READS", SSTR(par.numIterationsReads).c_str());

    //par.includeOnlyExtendable = false;
    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());


    // Read extension specific parameters
    par.includeOnlyExtendable = par.ancientIncludeOnlyExtendReads;
    par.kmerSize = par.ancientKmerSizeReads;
    par.kmersPerSequence = par.ancientkmersPerSequence;
    par.kmersPerSequenceScale = par.ancientKmersPerSequenceScale;
    cmd.addVariable("KMERMATCHER_READS_PAR", par.createParameterString(par.ancientKmermatcherReads).c_str());

    par.filterHits = false;
    // par.seqIdThr = par.corrContigSeqId;
    cmd.addVariable("UNGAPPED_ALN_PAR_ANCIENT_READS", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR_ANCIENT_READS", par.createParameterString(par.assembleresults).c_str());


    // Contig merging specific parameters
    par.includeOnlyExtendable = par.ancientIncludeOnlyExtendContigs;
    par.kmerSize = par.ancientKmerSizeContigs;
    par.kmersPerSequence = par.ancientkmersPerSequence;
    par.kmersPerSequenceScale = par.ancientKmersPerSequenceScale;
    par.seqIdThr = par.corrContigSeqId;
    cmd.addVariable("KMERMATCHER_CONTIGS_PAR", par.createParameterString(par.ancientKmermatcherContigs).c_str());
    cmd.addVariable("UNGAPPED_ALN_PAR_ANCIENT_CONTIGS", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR_ANCIENT_CONTIGS", par.createParameterString(par.assembleresults).c_str());


    cmd.addVariable("CALL_CYCLE_CHECK", par.cycleCheck ? "TRUE" : NULL);
    cmd.addVariable("CYCLE_CHECK_PAR", par.createParameterString(par.cyclecheck).c_str());

    cmd.addVariable("MIN_CONTIG_LEN", SSTR(par.minContigLen).c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("DB_MODE", par.dbMode ? "TRUE" : NULL);

    // Louis was here
    cmd.addVariable("RANDOM_ALN_PENAL", SSTR(par.randomAlignPenal).c_str());
    cmd.addVariable("PARAM_EXCESS_PENAL", SSTR(par.excessPenal).c_str());
    cmd.addVariable("CORRECTION_THRESHOLD", SSTR(par.corrReadsRySeqId).c_str());
    cmd.addVariable("CORRECTION_THRESHOLD_SEQID", SSTR(par.corrReadsSeqId).c_str());
    cmd.addVariable("LIKELIHOOD_THRESHOLD", SSTR(par.likelihoodThreshold).c_str());
    cmd.addVariable("SEQ_ID_MERGE_THRESH", SSTR(par.mergeSeqIdThr).c_str());
    cmd.addVariable("SEQ_ID_THRESH_CORR_CONTIG", SSTR(par.corrContigSeqId).c_str());
    cmd.addVariable("PATH_TO_DAMAGE_PATTERN", SSTR(par.ancientDamagePath).c_str());


    cmd.addVariable("KMER_SIZE_READS_ANCIENT", SSTR(par.ancientKmerSizeReads).c_str());
    cmd.addVariable("KMER_SIZE_CONTIGS_ANCIENT", SSTR(par.ancientKmerSizeContigs).c_str());
    cmd.addVariable("KMERS_PER_SEQ_ANCIENT", SSTR(par.ancientkmersPerSequence).c_str());
    cmd.addVariable("KMERS_PER_SEQ_SCALE_ANCIENT", SSTR(par.ancientKmersPerSequenceScale).c_str());
    cmd.addVariable("MIN_RYMER_SEQ_ID_ANCIENT", SSTR(par.rySeqIdThr).c_str());
    cmd.addVariable("UNSAFE_MODE", SSTR(par.ancientUnsafe).c_str());


    FileUtil::writeFile(tmpDir + "/nuclassemble.sh", nuclassemble_sh, nuclassemble_sh_len);
    std::string program(tmpDir + "/nuclassemble.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
