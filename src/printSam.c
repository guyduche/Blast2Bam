/*
The MIT License (MIT)

Copyright (c) 2015 Aurelien Guy-Duche

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "blastSam.h"

/************************************************************************************/
/*	Print SAM header																*/
/************************************************************************************/
static void sq_line(AppParamPtr app, char* filename) // Print the reference sequence dictionary (SQ lines)
{
	FILE* reader;
	char* str = NULL;
	int c, countSpace = 0;
	size_t lenStr = 0;

	reader = safeFOpen(filename, "r"); // Open the reference dictionary file (.dict)

	do
	{
		if (countSpace > 2 || !countSpace) // Skip the first line of the file and the end of the other lines
			do c = fgetc(reader); while (c != '\n' && c != EOF);

		lenStr = 0;
		countSpace = 0;
		str = NULL;
		c = fgetc(reader);
		lenStr++;

		while (c != '\n' && c != EOF && countSpace <= 2) // Keep only the reference name and its length
		{
			if (c == '\t') countSpace++;
			str = (char*) safeRealloc(str, lenStr+1);
			str[lenStr-1] = (char) c;
			c = fgetc(reader);
			lenStr++;
		}

		fwrite(str, sizeof(char), lenStr-1, app->out); // Print the SQ line in the SAM file
		free(str);
		if (c != EOF) fputc('\n', app->out);

	} while (c != EOF); // If there are more than one reference

	fclose(reader);
}

static int rg_line(AppParamPtr app, char* readGroup)
{	
	if (strstr(readGroup, "@RG") != readGroup) return 1; // Verify that the read group begins with @RG
	if (strstr(readGroup, "\tID:") == NULL) return 1; // Verify that the read group has an ID
	fprintf(app->out, "%s\n", readGroup); // Print the read group in the SAM header
	return 0;
}

int samHead(AppParamPtr app) // Print SAM header section
{
	sq_line(app, app->db); // Print the SQ Line
	if (app->readGroup != NULL) // If the read group is available and is correctly formatted, prints it
		if (rg_line(app,app->readGroup) == 1) return 1; // if the read group is incorrectly formatted, returns an error
	fprintf(app->out, "%s\n", app->pg_line); // Print the PG line
	return 0;
}


/************************************************************************************/
/*	Build the CIGAR string corresponding to the read alignment						*/
/************************************************************************************/
#define CSTRMACRO(TAG, FUNCT) { \
		COUNTCHARMACRO(FUNCT); \
		cigar->elements[cigar->size-1]->count = count; \
		cigar->elements[cigar->size-1]->symbol = TAG;}

#define COUNTCHARMACRO(FUNCT) { \
		count = 1; \
		pos++; \
		while (pos < (hsp->hsp_align_len) && FUNCT) \
			count++, pos++; \
		pos--;}

static CigarPtr cigarStrBuilding(SamOutputPtr samOut)
{
	int pos = 0, count = 0;
	size_t queryLength = samOut->query->read_len;
	HspPtr hsp = samOut->hsp;
	CigarPtr cigar = (CigarPtr) safeCalloc(1, sizeof(Cigar));

	if (hsp->hsp_query_from > 1) // 5' Soft clipping
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeCalloc(cigar->size, sizeof(CigarElementPtr)); // Create an array of elements for the CIGAR string
		cigar->elements[0] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement)); // Create the first element of the array
		cigar->elements[0]->count = hsp->hsp_query_from - 1; // Count = first alignment position in the query sequence - 1
		cigar->elements[0]->symbol = 'S';
	}

	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr)); // Create or append the array of elements of the CIGAR string
		cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement)); // Create a new element

		if (hsp->hsp_hseq[pos] == '-')
			CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-')) // Count the number of insertions

		else if (hsp->hsp_qseq[pos] == '-')
			CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-')) // Count the number of deletions

		else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
			CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])) // Count the number of matches

		else
			CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos])) // Count the number of mismatches

		if (cigar->elements[cigar->size-1]->count >= 100 && cigar->elements[cigar->size-1]->symbol == 'D')
			cigar->elements[cigar->size-1]->symbol = 'N'; // If there is more than a hundred deletion at a time, it is considered a skipped region

		if (cigar->elements[cigar->size-1]->symbol != '=')
			cigar->nbDiff += cigar->elements[cigar->size-1]->count; // Count the number of I, D/N and X for the NM tag
	}

	if ((queryLength - hsp->hsp_query_to) > 0) // 3' Soft clipping
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr));
		cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));
		cigar->elements[cigar->size-1]->count = queryLength - hsp->hsp_query_to; // Count = length of the read - position of the last aligned base on the query sequence 
		cigar->elements[cigar->size-1]->symbol = 'S';
	}

	return cigar;
}


/************************************************************************************/
/*	Get the reverse complement of a sequence										*/
/************************************************************************************/
static char* revStr(char* oldStr)
{
	size_t i = 0;
	int l = strlen(oldStr);
	char* newStr = safeCalloc(l+1, sizeof(char));

	for (--l; l >= 0; l--, i++)
	{
		switch (oldStr[l])
		{
			case 'A': case 'a' : newStr[i] = 'T'; break;
			case 'T': case 't' : newStr[i] = 'A'; break;
			case 'C': case 'c' : newStr[i] = 'G'; break;
			case 'G': case 'g' : newStr[i] = 'C'; break;
			default: newStr[i] = oldStr[l]; break;
		}
	}
	return newStr;
}


/************************************************************************************/
/*	Extract the position in the reference name, if any								*/
/************************************************************************************/
static int firstPosRef(const char* rname)
{
	char* colon = strchr(rname, ':');
	if (colon == NULL) 
		return 0;
	return strtol(colon + 1, NULL, 10);
}

/************************************************************************************/
/*	Print SAM alignment section														*/
/************************************************************************************/
#define SAM_PAIRED 0x1			// Paired end
#define SAM_PROPER_PAIR 0x2		// Read mapped in a proper pair
#define SAM_UNMAP 0x4			// Read unmapped
#define SAM_MUNMAP 0x8			// Mate unmapped
#define SAM_REVERSE 0x10		// Read mapped on the reverse strand
#define SAM_MREVERSE 0x20		// Mate mapped on the reverse strand
#define SAM_READF 0x40			// Read is first in pair
#define SAM_READL 0x80			// Read is last in pair
#define SAM_SECONDARY 0x100		// Not primary alignment
#define SAM_QCFAIL 0x200		// Failed control quality
#define SAM_DUP 0x400			// Optical or PCR duplicate
#define SAM_SUPPLEMENTARY 0x800	// Supplementary alignment

// Function to filter the results according to option -W
static int allowedToPrint(SamOutputPtr* samOut, int minLen, int countRec, int countHSPsec, int* countUnprint)
{
	if (samOut[0]->hsp != NULL && samOut[0]->hsp->hsp_align_len > minLen)
	{
		if (samOut[1] != NULL && samOut[1]->hsp != NULL && samOut[1]->hsp->hsp_align_len < minLen)
		{
			if (*countUnprint != (countRec - countHSPsec))
				{(*countUnprint)++; return 0;}	
			else
				samOut[1]->hsp = NULL;
		}
	}
	else
	{
		if (samOut[1] != NULL && samOut[1]->hsp != NULL && samOut[1]->hsp->hsp_align_len > minLen)
		{
			if (*countUnprint != (countRec - countHSPsec))
				{(*countUnprint)++; return 0;}
			else
				samOut[0]->hsp = NULL;
		}
		else
		{
			if (*countUnprint != countRec -1)
				{(*countUnprint)++; return 0;}
			else
			{
				samOut[0]->hsp = NULL;
				if (samOut[1] != NULL)
					samOut[1]->hsp = NULL;
			}
		}
	}
	return 1;
}

// Print the CIGAR str straight or reverse depending on the flag
static void printCigarStr (AppParamPtr app, CigarElementPtr* cigElements, size_t size, int flag)
{
	int i = 0;
	if (flag & SAM_REVERSE)
		for (i = size - 1; i >= 0; i--)
		{
			fprintf(app->out, "%d%c", cigElements[i]->count, cigElements[i]->symbol);
			free(cigElements[i]);
		}

	else
		for (i = 0; i < size; i++)
		{
			fprintf(app->out, "%d%c", cigElements[i]->count, cigElements[i]->symbol);
			free(cigElements[i]);
		}

	free(cigElements);
}

// Structure that contains the infos of one line of SAM alignment section
typedef struct SamLine
{
	char* readName;		// Col 1
	char* refName;		// Col 3
	char* rnext;		// Col 7
	char* seq;			// Col 10
	char* qual;			// Col 11
	int flag;			// Col 2
	int pos;			// Col 4
	int mapq;			// Col 5
	int pnext;			// Col 8
	int tlen;			// Col 9
	CigarPtr cigarStr;	// Col 6
} SamLine, *SamLinePtr;

// Print one line of the SAM alignment section
static void printSamLine (AppParamPtr app, SamLinePtr samLine)
{
	int i = 0;
	size_t qualLen = 0;
	char* seq = NULL;
	
	fputs(samLine->readName, app->out);
	fprintf(app->out, "\t%d\t", samLine->flag);
	fputs(samLine->refName, app->out);
	fprintf(app->out, "\t%d\t%d\t", samLine->pos, samLine->mapq);
	
	if (samLine->cigarStr != NULL)
	{
		printCigarStr (app, samLine->cigarStr->elements, samLine->cigarStr->size, samLine->flag);
		free(samLine->cigarStr);
	}
	else
		fputc('*', app->out);
		
	fputc('\t', app->out);
	fputs(samLine->rnext, app->out);
	fprintf(app->out, "\t%d\t%d\t", samLine->pnext, samLine->tlen);
	
	if (samLine->flag & SAM_REVERSE)
	{
		seq = revStr(samLine->seq);
		fputs(seq, app->out);
		fputc('\t', app->out);
		
		qualLen = strlen(samLine->qual);
		for (i = qualLen - 1; i >= 0; i--)
			fputc(samLine->qual[i], app->out);
	}
	else
	{
		fputs(samLine->seq, app->out);
		fputc('\t', app->out);
		fputs(samLine->qual, app->out);
	}
	
	if (samLine->cigarStr != NULL)
		fprintf(app->out, "\tNM:i:%d", samLine->cigarStr->nbDiff);

	if (app->readGroupID != NULL)
		fprintf(app->out, "\tRG:Z:%s", app->readGroupID);
				
	fputc('\n', app->out);
}

// Print the alignment section
void printSam(IterationSamPtr itSam, AppParamPtr app)
{
	int i = 0, j = 0, k = 0;
	int invk = 0, len0 = 0, len1 = 0, doNotPrint = 0, countUnprint = 0;
	SamOutputPtr samOut[2] = {NULL, NULL};
	SamLinePtr samLine = NULL;

	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				invk = (k ? 0 : 1);
				samOut[k] = itSam->samHits[i]->rsSam[j]->samOut[k];
				samOut[invk] = itSam->samHits[i]->rsSam[j]->samOut[invk];
				
				if (samOut[k] == NULL || doNotPrint) continue;

				if (app->minLen != 0 && !k)
					if (!allowedToPrint(samOut, app->minLen, itSam->samHits[i]->countRec, itSam->samHits[i]->countHSPsec, &countUnprint))
						{doNotPrint = 1; continue;}
				
				samLine = (SamLinePtr) safeCalloc(1, sizeof(SamLine));
				
				if (samOut[invk] != NULL) // Paired end
				{
					samLine->flag |= SAM_PAIRED | (!k ? SAM_READF : SAM_READL);
					
					if (samOut[invk]->hsp != NULL)
					{
						samLine->flag |= (samOut[invk]->hsp->hsp_hit_to < samOut[invk]->hsp->hsp_hit_from ? SAM_MREVERSE : 0);
						samLine->pnext = (samLine->flag & SAM_MREVERSE ? samOut[invk]->hsp->hsp_hit_to : samOut[invk]->hsp->hsp_hit_from);
						if (app->posOnChr) samLine->pnext += firstPosRef(samOut[invk]->rname);
						if (samOut[k]->hsp != NULL)
						{
							samLine->tlen = samOut[invk]->hsp->hsp_hit_from - samOut[k]->hsp->hsp_hit_from;
							if (samLine->tlen)
								samLine->tlen += (samLine->tlen > 0 ? 1 : -1);
							
							len0 = abs(samOut[k]->hsp->hsp_hit_to - samOut[k]->hsp->hsp_hit_from) + 1;
							len1 = abs(samOut[invk]->hsp->hsp_hit_to - samOut[invk]->hsp->hsp_hit_from) + 1;
							if (abs(samLine->tlen) > 3 * (len0 >= len1 ? len0 : len1) || abs(samLine->tlen) < (len0 >= len1 ? len1 : len0))
							{
								if (countUnprint != (itSam->samHits[i]->countRec -1))
									{doNotPrint = 1; countUnprint++; free(samLine); continue;}
								else
								{
									samLine->flag &= ~SAM_MREVERSE;
									samLine->flag |= SAM_MUNMAP;
									samLine->rnext = "*";
									samLine->pnext = 0;
									samOut[k]->hsp = NULL;
									samOut[invk]->hsp = NULL;
								}
							}
							else
							{
								samLine->rnext = "=";
								samLine->flag |= SAM_PROPER_PAIR;
							}
						}
						else
							samLine->rnext = shortName(samOut[invk]->rname);
					}
					else
					{
						samLine->flag |= SAM_MUNMAP;
						samLine->rnext = "*";
						samLine->pnext = 0;
					}
				}
				else // Single end
				{
					samLine->rnext = "*";
					samLine->pnext = 0;
				}

				samLine->readName = samOut[k]->query->name;
				samLine->seq = samOut[k]->query->seq;
				samLine->qual = samOut[k]->query->qual;
				
				if (samOut[k]->hsp != NULL)
				{
					if (j - countUnprint != 0)
						samLine->flag |= SAM_SECONDARY;
					samLine->flag |= (samOut[k]->hsp->hsp_hit_to < samOut[k]->hsp->hsp_hit_from ? SAM_REVERSE : 0);
					samLine->refName = shortName(samOut[k]->rname);
					samLine->pos = (samLine->flag & SAM_REVERSE ? samOut[k]->hsp->hsp_hit_to : samOut[k]->hsp->hsp_hit_from);
					if (app->posOnChr)
						samLine->pos += firstPosRef(samOut[k]->rname);
					samLine->cigarStr = cigarStrBuilding(samOut[k]);
					samLine->mapq = 60;
				}
				else
				{
					samLine->flag |= SAM_UNMAP;
					samLine->refName = "*";
				}
				
				printSamLine (app, samLine);
				free(samLine);
			}
		}
		countUnprint = 0;
	}
}

