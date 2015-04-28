
#ifndef SHORTREAD_H_INCLUDED
#define SHORTREAD_H_INCLUDED

#include <zlib.h>

typedef struct ShortRead
{
	char* name;
	char* seq;
	char* qual;
	size_t read_len;
} ShortRead, *ShortReadPtr;

ShortReadPtr shortReadNext(gzFile in);
void shortReadFree(ShortReadPtr ptr);

#endif // SHORTREAD_H_INCLUDED
