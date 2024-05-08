#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

//
extern int assemble(int argc, const char **argv, const Command &command);
extern int nuclassemble(int argc, const char** argv, const Command &command);
extern int guidedNuclAssemble(int argc, const char** argv, const Command &command);
extern int assembleresult(int argc, const char** argv, const Command &command);
extern int guidedassembleresults(int argc, const char** argv, const Command &command);
extern int ancientReadsResults(int argc, const char** argv, const Command &command);
extern int ancientContigsResults(int argc, const char** argv, const Command &command);
extern int correction(int argc, const char** argv, const Command &command);
extern int filternoncoding(int argc, const char** argv, const Command &command);
extern int mergereads(int argc, const char** argv, const Command &command);
extern int findassemblystart(int argc, const char** argv, const Command &command);
extern int cyclecheck(int argc, const char** argv, const Command &command);
extern int createhdb(int argc, const char** argv, const Command &command);
#endif
