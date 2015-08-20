/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum and Aurelien Guy-Duche

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/
#ifndef SHORTREAD_H_INCLUDED
#define SHORTREAD_H_INCLUDED

#include <zlib.h>

// Read infos extracted from the fastQ
typedef struct ShortRead
{
    char* name;         // Read name
    char* seq;          // Read sequence
    char* qual;         // Read quality
    size_t read_len;    // Read length
} ShortRead, *ShortReadPtr;


/************************************************************************************/
/*  Prototypes                                                                      */
/************************************************************************************/
gzFile initFastQ(int* fasta, char* filename);
ShortReadPtr shortReadNext(gzFile in, int fasta);
void shortReadFree(ShortReadPtr ptr);

#endif // SHORTREAD_H_INCLUDED
