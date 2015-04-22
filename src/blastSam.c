
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
	struct SamOutput_t* next;
} SamOutput, *SamOutputPtr;

int cigarStrBuilding(SamOutputPtr samOut, HspPtr hsp, int queryLength)
{
	int pos = 0;
	int count = 0;
	size_t sizeCStr = 0;
	int* cigarStr = samOut->cigar;
	
	if (queryLength > (hsp->hsp_align_len) && (hsp->hsp_query_from) > 1)
	{
		sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, sizeCStr, int, 1)
		cigarStr[0] = hsp->hsp_query_from - 1;
		cigarStr[1] = 'S';
	}
	
	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{	
		sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, sizeCStr, int, 1)
		
		if (hsp->hsp_hseq[pos] == '-')
			CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-'))
		
		else if (hsp->hsp_qseq[pos] == '-')
			CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-'))
			
		else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
			CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos]))
		
		else
			CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos]))
		
		if (cigarStr[sizeCStr-2] >= 100 && cigarStr[sizeCStr-1] == 'D')
			cigarStr[sizeCStr-1] = 'N';
		
	}
	
	if (queryLength > (hsp->hsp_align_len) && (queryLength - hsp->hsp_query_to) > 1)
	{
		sizeCStr += 2;
		STRBUILDINGMACRO(cigarStr, sizeCStr, int, 1)
		cigarStr[sizeCStr-2] = queryLength - hsp->hsp_query_to;
		cigarStr[sizeCStr-1] = 'S';
	}
	
	samOut->sizeCStr = sizeCStr;
	samOut->cigar = cigarStr;
	return 0;
}

void printRead(SamOutputPtr samOut)
{
	size_t i = 0;
	fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t", samOut->qname, samOut->flag, samOut->rname, samOut->pos, samOut->mapq);
	if (samOut->cigar != NULL)
		for (i = 0; i < samOut->sizeCStr; i++)
			fprintf(stdout,(i % 2 == 0 ? "%d":"%c"), samOut->cigar[i]);
	else
		fprintf(stdout, "*");
	fprintf(stdout, "\t%s\t%d\t%d\t%s\t%s\n", samOut->rnext, samOut->pnext, samOut->tlen, samOut->seq, samOut->qual);
}

static char* shortRefName(char* name)
{
	char* p = strpbrk(name," \t"); 
	if(p != 0)
		*p = 0;
	return name;
}

int printSAM(IterationPtr itFor, IterationPtr itRev, ShortReadPtr seqFor, ShortReadPtr seqRev)
{
	HitPtr hitCur = NULL;
	HspPtr hspFor = NULL;
	HspPtr hspRev = NULL;
	HspPtr hspBest = NULL;
	int posBestHsp = 0;
	SamOutputPtr samOutFor = NULL;
	SamOutputPtr samOutRev = NULL;
	int nbHSP = 0;
	int i = 0;
	
	samOutFor->qname = itFor->iteration_query_def;
	
	if (itRev == NULL && seqRev == NULL) // Single end
	{
		if (strcmp(seqFor->name, itFor->iteration_query_def))
			ERROR("Mismatch between fastq and blast output", 1)
		
		samOutFor->rnext = "*";
		samOutFor->pnext = 0;
		
		hitCur = itFor->iteration_hits;
		
		while (hitCur != NULL)
		{
			
			hspFor = hitCur->hit_hsps;
			hspCur = hspFor;
			
			if (hspFor == NULL)
				// samout specific for no hsp
			
			while (hspFor != NULL)
			{
				// allocation of samout relative to hsp
				nbHsp++;
				if (hspFor->hsp_score > hspBest->hsp_score)
				{
					hspBest = hspFor;
					posBestHsp = nbHsp-1;
				}
				
				hspFor = hspFor->next;
			}
			hitCur = hitCur->next;
		}

			while (hspCur != NULL)
			{
				if (hspCur->hsp_score > hspFor->hsp_score) // Print the highest scoring HSP
					hspFor = hspCur;
				hspCur = hspCur->next;
			}

			if (cigarStrBuilding(samOut, hspFor, itFor->iteration_query_len) != 0)
				ERROR("Error while building the cigar string", 1)
			
			if (i = 0)
			{
				samOutFor->flag = 1;
				samOutFor->seq = seqFor->seq;
				samOutFor->qual = seqFor->qual;
			}
			
			else
			{
				samOutFor->flag = 2;
				samOutFor->seq = "*";
				samOutFor->qual = "*";
			}
				
			
			printRead(itFor->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "*", 0, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);
			itFor->iteration_hits = itFor->iteration_hits->next;
			i++;
			}
		}
			
		else
			printRead(itFor->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqFor->seq, seqFor->qual);
		
		samOutFor->flag = 0; // TODO: Make the flag
		
		printRead(samOutFor);
	}
	
	else // Paired end
	{
		if (!strcmp(seqFor->name, itFor->iteration_query_def) && !strcmp(seqRev->name, itRev->iteration_query_def))
		{
			samOutRev = (SamOutputPtr) calloc(1, sizeof(SamOutput); // TODO: safe calloc
			
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

				if (cigarStrBuilding(samOut, hspFor, itFor->iteration_query_len) != 0)
					ERROR("Error while building the cigar string", 1)
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

				if (cigarStrBuilding(samOut, hspRev, itRev->iteration_query_len) != 0)
					ERROR("Error while building the cigar string", 1)
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
						if ((0 < abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to)) && (abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to) < 1000000)) // Find the most "compatible" Forward and Reverse HSP; TODO: Create a function to test the compatibility
						{
							if (cigarStrBuilding(samOut, hspFor, itFor->iteration_query_len) != 0)
								ERROR("Error while building the cigar string", 1)
							printRead(itFor->iteration_query_def, flag, shortRefName(itFor->iteration_hits->hit_def), hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", hspRev->hsp_hit_from, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);

							if (cigarStrBuilding(samOut, hspRev, itRev->iteration_query_len) != 0)
								ERROR("Error while building the cigar string", 1)
							printRead(itRev->iteration_query_def, flag, shortRefName(itRev->iteration_hits->hit_def), hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", hspFor->hsp_hit_from, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq, seqRev->qual);
						}

						hspRev = hspRev->next;
					}
					hspFor = hspFor->next;
				}
			}
		}
		printRead(samOutFor);
		printRead(samOutRev);
		
		if (samOutRev->cigar != NULL)
		free(samOutRev->cigar);
	free(samOutRev);
	}
	
	if (samOutFor->cigar != NULL)
		free(samOutFor->cigar);
	free(samOutFor);
	return 0;
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


ItSamPtr** hitRecord(HitPtr hitFor, HitPtr hitRev, ItSamPtr** itSam, RVarPtr rVar)
{	
	int i = 0;
	rVar->countRec++;
	
	itSam[rVar->countHit-1] = saferealloc(itSam[rVar->countHit-1], rVar->countRec * sizeof(ItSamPtr)); // Create and/or append a given Hit record table
	itSam[rVar->countHit-1][rVar->countRec-1] = safecalloc(1, sizeof(ItSam)); // Create a new record
	itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0] = safecalloc(1, sizeof(SamOutput)); // Create a structure to capture all the info concerning the forward strand
	itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0]->query = rVar->reads[0];
	
	if (hitRev == NULL) // Single End
	{	
		if (hitFor->hit_hsps == NULL) // Not mapped
		{
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0]->flag |= 0x4;
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0]->rname = "*";
			itSam[rVar->countHit-1][rVar->countRec-1]->score = 0;
		}
		else
		{
			// flag
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0]->rname = shortRefName(hitFor->hit_def);
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0]->hsp = hitFor->hit_hsps;
			if (cigarStrBuilding(itSam[rVar->countHit-1][rVar->countRec-1]->samOut[0], hitFor->hit_hsps, rVar->queryLength) != 0)
				ERROR("Error while building the cigar string", 1)
			itSam[rVar->countHit-1][rVar->countRec-1]->score = hitFor->hit_hsps->hsp_score;
			hitFor->hit_hsps = hitFor->hit_hsps->next;
			if (hitFor->hit_hsps != NULL)
				return hitRecord(hitFor, NULL, itSam, rVar); // Record the other HSPs if there are any
		}
		
		if (hitFor->next != NULL)
		{
			rVar->countHit++;
			rVar->countRec = 0;
			itSam = saferealloc(itSam, rVar->countHit * sizeof(ItSamPtr*)); // Append the Hit list
			return hitRecord(hitFor->next, NULL, itSam, rVar); // Record the other Hits if there are any
		}
		else
			return itSam;
	}
	
	// Only go there if it's paired end
	
	itSam[rVar->countHit-1][rVar->countRec-1]->samOut[1] = safecalloc(1, sizeof(SamOutput));
	itSam[rVar->countHit-1][rVar->countRec-1]->samOut[1]->query = rVar->reads[1];
	
	if (hitFor->hit_hsps == NULL && hitRev->hit_hsps == NULL) // Both unmapped
	{	
		for (i = 0; i < 2; i++)
		{
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[i]->flag |= 0x4;
			itSam[rVar->countHit-1][rVar->countRec-1]->samOut[i]->rname = "*";
		}
		
		itSam[rVar->countHit-1][rVar->countRec-1]->score = 0;
	}
	
	else if (

	

	// compare the name of the hits
	
		


// for the flag, to substract you can use flag &= ~0x100


typedef struct ItSam_t
{
	SamOutputPtr samOut[2];
	double score;
} ItSam, *ItSamPtr;

typedef struct RecordVariables_t
{
	int countHit;
	int countRec;
	ShortReadPtr reads[2];
} RVar, *RVarPtr;

typedef struct SamOutput_t
{
	ShortReadPtr query;
	HspPtr hsp;
	int flag;
	char* rname;
	int* cigar;
	size_t sizeCStr;
} SamOutput, *SamOutputPtr;


