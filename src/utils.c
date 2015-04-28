
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"

void* _safeMalloc(const char* file, int line, size_t n)
{
	void* p = malloc(n);
	if (p == NULL)
		fprintf(stderr, "Failed memory allocation in %s at line %d\n", file, line);
	return p;
}

void* _safeCalloc(const char* file, int line, size_t m, size_t n)
{
	void* p = calloc(m, n);
	if (p == NULL)
		fprintf(stderr, "Failed memory allocation in %s at line %d\n", file, line);
	return p;
}

void* _safeRealloc(const char* file, int line, void* ptr, size_t n)
{
	void* p = realloc(ptr, n);
	if (p == NULL)
		fprintf(stderr, "Failed memory allocation in %s at line %d\n", file, line);
	return p;
}

char* _safeStrdup(const char* file, int line, char* s)
{
	char* str = strdup(s);
	if (str == NULL)
		fprintf(stderr, "Failed memory allocation in %s at line %d\n", file, line);
	return str;
}

gzFile _safeGzOpen(const char* file, int line, char* filename, char* mode)
{
	gzFile fp = gzopen(filename, mode);
	if (fp == NULL)
		fprintf(stderr, "Unable to open the GZ file. Error in %s at line %d\n", file, line);
	return fp;
}

FILE* _safeFOpen(const char* file, int line, char* filename, char* mode)
{
	FILE* fp = fopen(filename, mode);
	if (fp == NULL)
		fprintf(stderr, "Unable to open the file. Error in %s at line %d\n", file, line);
	return fp;
}

xmlTextReaderPtr _safeXmlNewTextReaderFilename(const char* file, int line, char* filename)
{
	xmlTextReaderPtr fp = xmlNewTextReaderFilename(filename);
	if (fp == NULL)
		fprintf(stderr, "Unable to open the XML file. Error in %s at line %d\n", file, line);
	return fp;
}

int _safeXmlTextReaderRead(const char* file, int line, xmlTextReaderPtr fp)
{
	int evt = xmlTextReaderRead(fp);
	if (evt == -1)
		fprintf(stderr, "Error while reading the XML node. Error occurred in %s at line %d\n", file, line);
	return evt;
}






















