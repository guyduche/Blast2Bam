
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parseXML.h"
#include "blastSam.h"
#include "utils.h"
#include "shortRead.h"

#define CSTRMACRO(TAG, FUNCT) { \
		COUNTCHARMACRO(FUNCT); \
		cigarStr[*sizeCStr-2] = count; \
		cigarStr[*sizeCStr-1] = TAG;}
		
#define COUNTCHARMACRO(FUNCT) { \
		count = 1; \
		pos++; \
		while (pos < (hsp->hsp_align_len) && FUNCT) \
		{ \
			count++; \
			pos++; \
		} \
		pos--;}

typedef struct SamOutput_t
{
	char* qname;
	int flag;
	char* rname;
	int pos;
	int mapq;
	int* cigar;
	size_t sizeCStr;
	char* rnext;
	int pnext;
	int tlen;
	char* seq;
	char* qual;
} SamOutput, *SamOutputPtr;

int* cigarStrBuilding(int* cigarStr, Hsp* hsp, int queryLength, size_t* sizeCStr)
{
	int pos = 0;
	int count = 0;
	*sizeCStr = 0;
	
	if (queryLength > (hsp->hsp_align_len) && (hsp->hsp_query_from) > 1)
	{
		*sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		cigarStr[0] = hsp->hsp_query_from - 1;
		cigarStr[1] = 'S';
	}
	
	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{	
		*sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		
		if (hsp->hsp_hseq[pos] == '-')
			CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-'))
		
		else if (hsp->hsp_qseq[pos] == '-')
			CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-'))
			
		else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
			CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos]))
		
		else
			CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos]))
		
		if (cigarStr[*sizeCStr-2] >= 100 && cigarStr[*sizeCStr-1] == 'D')
			cigarStr[*sizeCStr-1] = 'N';
		
	}
	
	if (queryLength > (hsp->hsp_align_len) && (queryLength - hsp->hsp_query_to) > 1)
	{
		*sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		cigarStr[*sizeCStr-2] = queryLength - hsp->hsp_query_to;
		cigarStr[*sizeCStr-1] = 'S';
	}
	
	return cigarStr;
}

void printRead(char* qname, int flag, char* rname, int pos, int mapq, int* cigar, size_t sizeCStr, char* rnext, int pnext, int tlen, char* seq, char* qual)
{
	size_t i = 0;
	fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t", qname, flag, rname, pos, mapq);
	if (cigar != NULL)
		for (i = 0; i < sizeCStr; i++)
			fprintf(stdout,(i % 2 == 0 ? "%d":"%c"), cigar[i]);
	else
		fprintf(stdout, "*");
	fprintf(stdout, "\t%s\t%d\t%d\t%s\t%s\n", rnext, pnext, tlen, seq, qual);
}

static char* shortRefName(char* name)
{
	char* p = strpbrk(name," \t"); 
	if(p != 0)
		*p = 0;
	return name;
}

void printSAM(IterationPtr itFor, IterationPtr itRev, ShortReadPtr seqFor, ShortReadPtr seqRev)
{
	HspPtr hspFor = NULL;
	HspPtr hspRev = NULL;
	HspPtr hspCur = NULL;
	int* cigar = NULL;
	size_t sizeCStr = 0;
	int flag = 0;
	
	if (itRev == NULL && seqRev == NULL) // Single end
	{
		if (!strcmp(seqFor->name, itFor->iteration_query_def))
		{
			if (itFor->iteration_hits->hit_hsps != NULL)
			{
				hspFor = itFor->iteration_hits->hit_hsps;
				hspCur = hspFor;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspFor->hsp_score) // Print the highest scoring HSP
						hspFor = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hspFor, itFor->iteration_query_len, &sizeCStr);
				printRead(itFor->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "*", 0, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);
			}
			
			else
				printRead(itFor->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqFor->seq, seqFor->qual);
		}
	}
	
	else // Paired end
	{
		if (!strcmp(seqFor->name, itFor->iteration_query_def) && !strcmp(seqRev->name, itRev->iteration_query_def))
		{
			if (itFor->iteration_hits->hit_hsps == NULL && itRev->iteration_hits->hit_hsps == NULL) // Not mapped
			{
				printRead(itFor->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqFor->seq, seqFor->qual);
				printRead(itRev->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqRev->seq, seqRev->qual);
			}

			else if (itRev->iteration_hits->hit_hsps == NULL) // Forward is mapped, Reverse is not
			{
				hspFor = itFor->iteration_hits->hit_hsps;
				hspCur = hspFor;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspFor->hsp_score) // Print the highest scoring Forward HSP
						hspFor = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hspFor, itFor->iteration_query_len, &sizeCStr);
				printRead(itFor->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", 0, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);

				printRead(itRev->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), 0, 0, NULL, 0, "=", hspFor->hsp_hit_from, 0, seqRev->seq, seqRev->qual);
			}

			else if (itFor->iteration_hits->hit_hsps == NULL) // Reverse is mapped, Forward is not
			{
				hspRev = itRev->iteration_hits->hit_hsps;
				hspCur = hspRev;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspRev->hsp_score) // Print the highest scoring Reverse HSP
						hspRev = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hspRev, itRev->iteration_query_len, &sizeCStr);
				printRead(itRev->iteration_query_def, flag, shortRefName(itRev->iteration_hits->hit_def), hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", 0, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq, seqRev->qual);

				printRead(itFor->iteration_query_def, flag, shortRefName(itRev->iteration_hits->hit_def), 0, 0, NULL, 0, "=", hspRev->hsp_hit_from, 0, seqFor->seq, seqFor->qual);
			}

			else // both are mapped
			{			
				hspFor = itFor->iteration_hits->hit_hsps;
				hspRev = itRev->iteration_hits->hit_hsps;
				
				while (hspFor != NULL)
				{
					hspRev = itRev->iteration_hits->hit_hsps;
					while (hspRev != NULL)
					{
						if ((0 < abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to)) && (abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to) < 1000000)) // Find the most "compatible" Forward and Reverse HSP; NEED SOME WORK !!!
						{
							cigar = cigarStrBuilding(cigar, hspFor, itFor->iteration_query_len, &sizeCStr);
							printRead(itFor->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", hspRev->hsp_hit_from, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);

							cigar = cigarStrBuilding(cigar, hspRev, itRev->iteration_query_len, &sizeCStr);
							printRead(itRev->iteration_query_def, flag, shortRefName(itRev->iteration_hits->hit_def), hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", hspFor->hsp_hit_from, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq, seqRev->qual);
						}

						hspRev = hspRev->next;
					}
					hspFor = hspFor->next;
				}
			}
		}
	}
	
	if (cigar != NULL)
		free(cigar);
}

int parseDict(char* filename)
{
	FILE* reader;
	char* str = NULL;
	int c;
	int lenStr = 0;
	int countSpace = 0;
	
	reader = fopen(filename, "r");

	if (reader == NULL)
		ERROR("Error while opening the file\n", 1)

	do
	{
		if (countSpace > 2 || !countSpace)
			do
			{
				c = fgetc(reader);
			} while (c != '\n' && c != EOF);

		lenStr = 0;
		countSpace = 0;
		str = NULL;
		c = fgetc(reader);
		lenStr++;

		while (c != '\n' && c != EOF && countSpace <= 2)
		{
			if (c == '\t')
				countSpace++;			
			
			STRBUILDINGMACRO(str, lenStr+1, char, 1)
			str[lenStr-1] = (char) c;
			c = fgetc(reader);
			lenStr++;
		}

		STRBUILDINGMACRO(str, lenStr, char, 1)
		str[lenStr-1] = '\0';

		fprintf(stdout, "%s", str);
		free(str);
		if (c != EOF)
			fprintf(stdout, "\n");

	} while (c != EOF);
	
	fclose(reader);
	
	return 0;
}




