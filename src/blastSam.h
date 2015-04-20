
#ifndef SAM_H_INCLUDED
#define SAM_H_INCLUDED

#include "parseXML.h"
#include "shortRead.h"

int* cigarStrBuilding(int* cigarStr, Hsp* hsp, int queryLength, size_t* sizeCStr);
void printRead(char* qname, int flag, char* rname, int pos, int mapq, int* cigar, size_t sizeCStr, char* rnext, int pnext, int tlen, char* seq, char* qual);
char* shortRefName(char* name);
void printSAM(IterationPtr itFor, IterationPtr itRev, ShortReadPtr seqFor, ShortReadPtr seqRev);
int parseDict(char* filename);

#endif // SAM_H_INCLUDED
