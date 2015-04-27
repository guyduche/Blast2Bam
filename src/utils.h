
#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <stdio.h>
#include <libxml/xmlreader.h>
#include <zlib.h>

#define ERROR(MESSAGE, ENDING) { \
		fprintf(stderr, MESSAGE); \
		return ENDING;}


#define safeMalloc(n) _safeMalloc(__FILE__, __LINE__, n)
#define safeCalloc(m, n) _safeCalloc(__FILE__, __LINE__, m, n)
#define safeRealloc(ptr, n) _safeRealloc(__FILE__,__LINE__, ptr, n)
#define safeStrdup(s) _safeStrdup(__FILE__, __LINE__, s)
#define safeGzOpen(filename, mode) _safeGzOpen(__FILE__, __LINE__, filename, mode)
#define safeFOpen(filename, mode) _safeFOpen(__FILE__, __LINE__, filename, mode)
#define safeXmlNewTextReaderFilename(filename) _safeXmlNewTextReaderFilename(__FILE__, __LINE__, filename)
#define safeXmlTextReaderRead(fp) _safeXmlTextReaderRead (__FILE__, __LINE__, fp)


void* _safeMalloc(const char* file, int line, size_t n);
void* _safeCalloc(const char* file, int line, size_t m, size_t n);
void* _safeRealloc(const char* file, int line, void* ptr, size_t n);
char* _safeStrdup(const char* file, int line, char* s);
gzFile _safeGzOpen(const char* file, int line, char* filename, char* mode);
FILE* _safeFOpen(const char* file, int line, char* filename, char* mode);
xmlTextReaderPtr _safeXmlNewTextReaderFilename(const char* file, int line, char* filename);
int _safeXmlTextReaderRead(const char* file, int line, xmlTextReaderPtr fp);

#endif // UTILS_H_INCLUDED
