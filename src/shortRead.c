#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "utils.h"
#include "shortRead.h"


typedef struct ShortRead
{
	char* name;
	char* seq;
	char* qual;
	size_t read_len;
} ShortRead, *_ShortReadPtr;

#define CAST(a) ((_ShortReadPtr)(a))

static char* gzReadLine(gzFile in, size_t* line_len)
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
	
	line = (char*) safeMalloc(len);
	memcpy(line, buffer, len);
	
	return line;
}

void shortReadFree(ShortReadPtr ptr)
{
	if (ptr == NULL) return;
	free(CAST(ptr)->name - 1);
	free(CAST(ptr)->seq);
	free(CAST(ptr)->qual);
	free(ptr);
}

ShortReadPtr shortReadNext(gzFile in)
{
	size_t line_len = 0UL;
	ShortReadPtr ptr = NULL;
	
	if (gzeof(in)) return NULL;
	
	ptr = (ShortReadPtr) safeCalloc(1, sizeof(ShortRead));
	
	ptr->name = gzReadLine(in, &line_len) + 1;
	ptr->seq = gzReadLine(in, &line_len);
	ptr->read_len = line_len;
	free(gzReadLine(in, &line_len)); // Discard the 3rd line
	ptr->qual = gzReadLine(in, &line_len);
	
	if (line_len != ptr->read_len)
		ERROR("Wrong quality string length", NULL)
		
	return ptr;
}

