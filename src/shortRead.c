#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "macro.h"
#include "shortRead.h"

ShortReadPtr ShortReadNew()
{
	ShortReadPtr p = (ShortReadPtr) calloc(1, sizeof(ShortRead));
	if (p == NULL)
		ERROR("Failed memory allocation of the new short read", NULL)
	return p;
}

void ShortReadFree(ShortReadPtr ptr)
{
	if (ptr == NULL) return;
	free(ptr->name - 1);
	free(ptr->seq);
	free(ptr->qual);
	free(ptr);
}

char* gzReadLine(gzFile in, size_t* line_len)
{
	size_t len;
	char buffer[BUFSIZ];
	char* line = NULL;
	
	if (gzeof(in)) return NULL;
	if (gzgets(in, buffer, BUFSIZ) == NULL) return NULL; // Warning verify end of string
	
	len = strlen(buffer);
	if (len == 0 || buffer[len-1] != '\n')
		ERROR("Error while reading the new line", NULL)
	
	buffer[len-1] = '\0';
	*line_len = len-1;
	
	line = (char*) malloc(len);
	if (line == NULL)
		ERROR("Failed memory allocation of the new line", NULL)
	memcpy(line, buffer, len);
	
	return line;
}

ShortReadPtr ShortReadNext(gzFile in)
{
	size_t line_len = 0UL;
	ShortReadPtr ptr = NULL;
	
	if (gzeof(in)) return NULL;
	
	ptr = ShortReadNew();
	
	ptr->name = gzReadLine(in, &line_len) + 1;
	ptr->seq = gzReadLine(in, &line_len);
	ptr->read_len = line_len;
	free(gzReadLine(in, &line_len)); // Discard the 3rd line
	ptr->qual = gzReadLine(in, &line_len);
	
	if (line_len != ptr->read_len)
		ERROR("Wrong quality string length", NULL)
		
	return ptr;
}

