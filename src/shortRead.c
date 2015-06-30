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
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "shortRead.h"

/************************************************************************************/
/*  FastQ parsing                                                                   */
/************************************************************************************/
// Get one line of the file
static char* gzReadLine(gzFile in, size_t* line_len)
{
    size_t len = 0;
    int endsWithCR = 0;
    char* buffer = NULL, *line = NULL;
    *line_len = 0;

    while (!endsWithCR) // Until a '\n' is found
    { // Useful in case BUFSIZ is too small to get the whole line
        buffer = (char*) safeMalloc(BUFSIZ);
        if (gzgets(in, buffer, BUFSIZ) == NULL) return NULL; 
        if (gzeof(in)) return NULL;

        len = strlen(buffer);
        if (!len) ERROR("Error while reading the new line\n", NULL);
        if (buffer[len-1] == '\n')
        {
            buffer[len-1] = '\0';
            endsWithCR = 1;
        }

        line = (char*) safeRealloc(line, *line_len + len + 1);
        if (*line_len == 0) *line = '\0';   // Initialize line to enable the use of strncat
        strncat(line, buffer, len);

        *line_len += len;
        free(buffer);
    }
    (*line_len)--;
    return line;
}

// Get the infos of one read
ShortReadPtr shortReadNext(gzFile in)
{
    size_t line_len = 0;
    ShortReadPtr ptr = NULL;

    ptr = (ShortReadPtr) safeCalloc(1, sizeof(ShortRead));

    ptr->name = gzReadLine(in, &line_len);          // Get the read name
    if (ptr->name == NULL)
        return NULL;
    else
        ptr->name = shortName(ptr->name + 1);       // Put the short version of the read name in the ShortRead structure
    ptr->seq = gzReadLine(in, &line_len);           // Get the sequence
    ptr->read_len = line_len;

    switch (gzgetc(in))
    {
        case -1: break;                             // End of file
        
        case '+':                                   // FastQ
        {
            free(gzReadLine(in, &line_len));            // Discard the third line
            ptr->qual = gzReadLine(in, &line_len);      // Get the sequence quality
            if (line_len != ptr->read_len)              // Seq and qual should have the same length
                ERROR("Wrong quality string length\n", NULL);
        } break;
        
        case '>':                                   // Fasta
            gzseek(in, -1, SEEK_CUR); break;
        
        default:
            ERROR("Neither FastQ nor Fasta\n", NULL); break;
    }

    return ptr;
}

/************************************************************************************/
/*  Deallocation of the shortRead structure                                         */
/************************************************************************************/
void shortReadFree(ShortReadPtr ptr)
{
    if (ptr == NULL) return;
    free(ptr->name - 1);
    free(ptr->seq);
    free(ptr->qual);
    free(ptr);
}
