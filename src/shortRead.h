
#ifndef SHORTREAD_H_INCLUDED
#define SHORTREAD_H_INCLUDED

#include "parseXML.h"

typedef struct shortRead_t
{
	char* name;
	char* seq;
	char* qual;
	size_t read_len;
} ShortRead, *ShortReadPtr;

typedef struct ShortReadAndIteration
{
	ShortReadPtr seq;
	IterationPtr blast;
} ShortReadAndIteration;

ShortReadPtr ShortReadNew();
char* gzReadLine(gzFile in, size_t* line_len);
ShortReadPtr ShortReadNext(gzFile in);
void ShortReadFree(ShortReadPtr ptr);

#endif // SHORTREAD_H_INCLUDED
