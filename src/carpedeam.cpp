#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "carpedeam";
const char* tool_name = "CarpeDeam";
const char* tool_introduction = "Ancient Metagenome Assembler.";
const char* main_author = "Louis Kraft (loipwr3000@gmail.com)";
const char* show_extended_help = NULL;
const char* show_bash_info = NULL;
bool hide_base_commands = true;
void (*validatorUpdate)(void) = 0;
LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        // Plass workflows
        {"ancient_assemble",         guidedNuclAssemble,            &localPar.guidedNuclAssembleworkflow, COMMAND_MAIN,
                "Damage aware nucleotide assembly iterative greedy overlap assembly using damage matrix and nucleotide information => CarpeDeam",
                NULL,
                "Louis Kraft <lokraf@dtu.dk>",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                "Kraft L et al.: A De Novo Metagenome Assembler for Heavily Damaged Ancient Datasets (published soon)", {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"nuclassemble",          nuclassemble,     &localPar.nuclassembleworkflow, COMMAND_MAIN,
                "Modified nuclassemble module from PenguiN",
                NULL,
                "Louis Kraft <lokraf@dtu.dk> and Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        // Plass tools
        {"ancient_read_assemble",      ancientReadsResults,       &localPar.assembleresults,      COMMAND_HIDDEN,
                "Extending representative sequence to the left and right side using ungapped alignments.",
                NULL,
                "Louis Kraft <lokraf@dtu.dk>",
                "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                        {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                        {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"ancient_contig_merge",      ancientContigsResults,       &localPar.assembleresults,      COMMAND_HIDDEN,
                "Extending representative sequence to the left and right side using ungapped alignments.",
                NULL,
                "Louis Kraft <lokraf@dtu.dk>",
                "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                        {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                        {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"ancient_correction",      correction,       &localPar.assembleresults, COMMAND_MAIN,
                "Correct deaminations",
                NULL,
                "Louis Kraft <lokraf@dtu.dk>",
                "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                        {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                        {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"mergereads",      mergereads,      &localPar.onlythreads,          COMMAND_HIDDEN,
                "Merge paired-end reads from FASTQ file (powered by FLASH)",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <o:sequenceDB>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"cyclecheck",      cyclecheck,      &localPar.cyclecheck,          COMMAND_HIDDEN,
                "Simple cycle detector",
                NULL,
                "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
                "<i:sequenceDB> <o:sequenceDBcycle>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                {"cycleResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::nuclDb }}},
        {"createhdb",      createhdb,      &localPar.createhdb,          COMMAND_HIDDEN,
                "Generate header db file for given sequence db file",
                NULL,
                "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
                "<i:sequenceDB> [<i:sequenceDBcycle>] <o:headerDB>",
                CITATION_PLASS, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};
