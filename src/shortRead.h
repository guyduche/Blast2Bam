
#ifndef SHORTREAD_H_INCLUDED
#define SHORTREAD_H_INCLUDED

#include <zlib.h>
#include "parseXML.h"

typedef void* ShortReadPtr;

ShortReadPtr shortReadNext(gzFile in);
void shortReadFree(ShortReadPtr ptr);
void const char* shortReadName(ShortReadPtr);
void const char* shortReadSequence(ShortReadPtr);
void const char* shortReadQuality(ShortReadPtr);
void int shortReadSize(ShortReadPtr);

#endif // SHORTREAD_H_INCLUDED
