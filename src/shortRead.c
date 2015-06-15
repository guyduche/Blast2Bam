#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "shortRead.h"

static char* gzReadLine(gzFile in, size_t* line_len)
{
	size_t len = 0;
	int endsWithCR = 0;
	char* buffer = NULL;
	char* line = NULL;
	*line_len = 0;

	while (!endsWithCR)
	{
		buffer = (char*) safeMalloc(BUFSIZ);
		if (gzgets(in, buffer, BUFSIZ) == NULL) return NULL; // Warning verify end of string
		if (gzeof(in)) return NULL;

		len = strlen(buffer);
		if (!len) ERROR("Error while reading the new line\n", NULL);
		if (buffer[len-1] == '\n')
		{
			buffer[len-1] = '\0';
			endsWithCR = 1;
		}

		line = (char*) safeRealloc(line, *line_len + len + 1);
		if (*line_len == 0) *line = '\0';
		strncat(line, buffer, len);

		*line_len += len;
		free(buffer);
	}
	(*line_len)--;
	return line;
}

void shortReadFree(ShortReadPtr ptr)
{
	if (ptr == NULL) return;
	free(ptr->name - 1);
	free(ptr->seq);
	free(ptr->qual);
	free(ptr);
}

ShortReadPtr shortReadNext(gzFile in)
{
	size_t line_len = 0;
	ShortReadPtr ptr = NULL;

	ptr = (ShortReadPtr) safeCalloc(1, sizeof(ShortRead));

	ptr->name = gzReadLine(in, &line_len);
	if (ptr->name == NULL)
		return NULL;
	else
		ptr->name = shortName(ptr->name + 1);
	ptr->seq = gzReadLine(in, &line_len);
	ptr->read_len = line_len;
	free(gzReadLine(in, &line_len)); // Discard the 3rd line
	ptr->qual = gzReadLine(in, &line_len);

	if (line_len != ptr->read_len)
		ERROR("Wrong quality string length\n", NULL);

	return ptr;
}

