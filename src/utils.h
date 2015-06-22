/*
The MIT License (MIT)

Copyright (c) 2015 Aurelien Guy-Duche and Pierre Lindenbaum

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
#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <stdio.h>
#include <libxml/xmlreader.h>
#include <zlib.h>


/************************************************************************************/
/*  Error macros                                                                    */
/************************************************************************************/
#define DEBUG(FORMAT, ARGS...) \
    do { \
    fprintf(stderr, "[%s:%d]:", __FILE__, __LINE__); \
    fprintf(stderr, FORMAT, ARGS); \
    fputc('\n', stderr); \
    } while(0)

#define ERROR(MESSAGE, ENDING) \
    do { \
    DEBUG("%s", MESSAGE); \
    return ENDING; \
    } while(0)


/************************************************************************************/
/*  Safe memory allocation macros                                                   */
/************************************************************************************/
#define safeMalloc(size) _safeMalloc(__FILE__, __LINE__, size)
#define safeCalloc(number, size) _safeCalloc(__FILE__, __LINE__, number, size)
#define safeRealloc(ptr, size) _safeRealloc(__FILE__,__LINE__, ptr, size)
#define safeStrdup(s) _safeStrdup(__FILE__, __LINE__, s)
#define safeStrAppend(x,y) _safeStrAppend(__FILE__, __LINE__, x, y)
#define safeGzOpen(filename, mode) _safeGzOpen(__FILE__, __LINE__, filename, mode)
#define safeFOpen(filename, mode) _safeFOpen(__FILE__, __LINE__, filename, mode)
#define safeXmlNewTextReaderFilename(filename) _safeXmlNewTextReaderFilename(__FILE__, __LINE__, filename)
#define safeXmlTextReaderRead(fp) _safeXmlTextReaderRead (__FILE__, __LINE__, fp)


/************************************************************************************/
/*  Prototypes                                                                      */
/************************************************************************************/
void* _safeMalloc(const char* file, int line, size_t size);
void* _safeCalloc(const char* file, int line, size_t number, size_t size);
void* _safeRealloc(const char* file, int line, void* ptr, size_t size);
char* _safeStrdup(const char* file, int line, char* s);
char* _safeStrAppend(const char* file, int line, char* x, const char* y);
gzFile _safeGzOpen(const char* file, int line, char* filename, char* mode);
FILE* _safeFOpen(const char* file, int line, char* filename, char* mode);
xmlTextReaderPtr _safeXmlNewTextReaderFilename(const char* file, int line, char* filename);
int _safeXmlTextReaderRead(const char* file, int line, xmlTextReaderPtr fp);
char* shortName(char* name);

#endif // UTILS_H_INCLUDED
