
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

/*
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
*/

typedef struct SamOutput_t
{
	ShortReadPtr query;
	HspPtr hsp;
	int flag;
	char* rname;
	int* cigar;
	size_t sizeCStr;
} SamOutput, *SamOutputPtr;

typedef struct RecordSam_t
{
	SamOutputPtr samOut[2];
	double score;
} RecordSam, *RecordSamPtr;

typedef struct SamHit_t
{
	RecordSamPtr* rsSam;
	size_t countRec;
} SamHit, *SamHitPtr;

typedef struct IterationSam_t
{
	SamHitPtr* samHits;
	size_t countHit;
} IterationSam, *IterationSamPtr;

typedef struct RecordVariables_t
{
	int end;
	int tmpHitNb;
	HitPtr hitTmp;
	HspPtr hspTmp;
	ShortReadPtr reads[2];
} RVar, *RVarPtr;



int cigarStrBuilding(SamOutputPtr samOut)
{
	int pos = 0;
	int count = 0;
	size_t sizeCStr = 0;
	int* cigarStr = samOut->cigar;
	HspPtr hsp = samOut->hsp;
	int queryLength = samOut->query->read_len;
	
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

void printSam(IterationSamPtr itSam)
{
	int i = 0;
	int j = 0;
	int k = 0;
	size_t l = 0;
	RecordSamPtr rSam;
	char* rnext = NULL;
	int pnext = 0;
	
	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				rSam = itSam->samHits[i]->rsSam[j];
				if (rSam->samOut[k] == NULL) continue;
				
				if (!k ? rSam->samOut[1] : rSam->samOut[0]) != NULL)
				{
					rnext = ((!k ? rSam->samOut[1]->rname : rSam->samOut[0]->rname) == NULL ? "*" : (!k ? rSam->samOut[1]->rname : rSam->samOut[0]->rname));
					pnext = ((!k ? rSam->samOut[1]->hsp : rSam->samOut[0]->hsp) == NULL ? 0 : (!k ? rSam->samOut[1]->hsp->hsp_hit_from : rSam->samOut[0]->hsp->hsp_hit_from));
				}
				else
				{
					rnext = "*";
					pnext = 0;
				}
					
				if (rSam->samOut[k]->hsp == NULL)
					fprintf(stdout, "%s\t%d\t*\t0\t0\t*\t%s\t%d\t0\t%s\t%s\n", rSam->samOut[k]->query->name, rSam->samOut[k]->flag, rnext, pnext, rSam->samOut[k]->query->seq, rSam->samOut[k]->query->qual);
				else
				{
					fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t", rSam->samOut[k]->query->name, rSam->samOut[k]->flag, rSam->samOut[k]->rname, rSam->samOut[k]->hsp->hsp_hit_from, rSam->samOut[k]->hsp->hsp_score);
					for (l = 0; l < rSam->samOut[k]->sizeCStr; l++)
						fprintf(stdout,(l % 2 == 0 ? "%d":"%c"), samOut->cigar[l]);
					fprintf(stdout, "\t%s\t%d\t%d\t%s\t%s\n", rnext, pnext, rSam->samOut[k]->hsp->hsp_hit_to - rSam->samOut[k]->hsp->hsp_hit_from, rSam->samOut[k]->query->seq, rSam->samOut[k]->query->qual);
				}
			}
		}
	}
}

static char* shortRefName(char* name)
{
	char* p = strpbrk(name," \t"); 
	if(p != 0)
		*p = 0;
	return name;
}

IterationSamPtr hitRecord(HitPtr hitFor, HitPtr hitRev, IterationSamPtr itSam, RVarPtr rVar)
{	
	int i = 0;
	size_t countHit = itSam->countHit;
	size_t countRec;
	RecordSamPtr rSam;
	
	itSam->samHits[countHit-1]->countRec++;
	countRec = itSam->samHits[countHit-1]->countRec;
	
	itSam->samHits[countHit-1]->rsSam = (RecordSamPtr*) saferealloc(itSam->samHits[countHit-1]->rsSam, countRec * sizeof(RecordSamPtr)); // Create and/or append a given Hit record table
	itSam->samHits[countHit-1]->rsSam[countRec-1] = (RecordSamPtr) safecalloc(1, sizeof(RecordSam)); // Create a new record
	rSam = itSam->samHits[countHit-1]->rsSam[countRec-1];
	
	if (!rVar->end)
	{
		rSam->samOut[0] = (SamOutputPtr) safecalloc(1, sizeof(SamOutput)); // Create a structure to capture all the info concerning the forward strand
		rSam->samOut[0]->query = rVar->reads[0];
	}
	
	if (hitRev != NULL)
	{
		rSam->samOut[1] = (SamOutputPtr) safecalloc(1, sizeof(SamOutput));
		rSam->samOut[1]->query = rVar->reads[1];
	}
	
	if (hitRev == NULL) // Single end or the forward strand is mapped on a reference the reverse strand is not
	{	
		if (hitFor->hit_hsps == NULL) // Not mapped
		{
			rSam->samOut[0]->flag |= 0x4;
			rSam->score = 0.0;
		}
		else
		{
			// TODO: flag
			rSam->samOut[0]->rname = shortRefName(hitFor->hit_def);
			rSam->samOut[0]->hsp = hitFor->hit_hsps;
			if (cigarStrBuilding(rSam->samOut[0]) != 0)
				ERROR("Error while building the cigar string", NULL)
			rSam->score = (double) hitFor->hit_hsps->hsp_score;
			hitFor->hit_hsps = hitFor->hit_hsps->next;
			if (hitFor->hit_hsps != NULL)
				return hitRecord(hitFor, NULL, itSam, rVar); // Record the other HSPs if there are any
		}
	}
	
	else if (rVar->end) // The reverse strand is mapped on a reference where the forward strand is not
	{
		if (rVar->tmpHitNb == countHit)
		{
			// TODO: flag
			rSam->samOut[1]->rname = shortRefName(hitRev->hit_def);
			rSam->samOut[1]->hsp = hitRev->hit_hsps;
			if (cigarStrBuilding(rSam->samOut[1]) != 0)
				ERROR("Error while building the cigar string", NULL)
			rSam->score = hitRev->hit_hsps->hsp_score;
			hitRev->hit_hsps = hitRev->hit_hsps->next;
			if (hitRev->hit_hsps != NULL)
				return hitRecord(hitFor, hitRev, itSam, rVar); // Record the other HSPs if there are any
		}
	}
	
	else if (hitFor->hit_def == NULL && hitRev->hit_def == NULL) // Both unmapped
	{	
		for (i = 0; i < 2; i++)
			rSam->samOut[i]->flag |= 0x4;
		
		rSam->score = 0.0;
	}
	
	else if (strcmp(hitFor->hit_def, hitRev->hit_def))
	{
		if (hitRev->next == NULL || !strcmp(hitFor->hit_def, hitRev->next->hit_def))
		{
			itSam->countHit++;
			itSam->samHits = (SamHitPtr*) saferealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
			itSam->samHits[itSam->countHit-1]->countRec = 0;
		}
		
		return hitRecord(hitFor, hitRev->next, itSam, rVar);
	}
	
	else if (strcmp(hitFor->hit_def, hitRev->hit_def) == 0) // Forward and reverse strand are mapped on the same reference
	{
		if (rVar->hspTmp == NULL)
			hspTmp = hitRev->hit_hsps;
		
		for (i = 0; i < 2; i++)
		{
			rSam->samOut[i]->hsp = (!i ? hitFor->hit_hsps : hitRev->hit_hsps);
			rSam->samOut[i]->rname = (!i ? hitFor->hit_def : hitRev->hit_def);
			if (cigarStrBuilding(rSam->samOut[i]) != 0)
				ERROR("Error while building the cigar string", NULL)
		}
		// TODO: flag and score
		hitRev->hit_hsps = hitRev->hit_hsps->next;
		if (hitRev->hit_hsps != NULL)
			return hitRecord(hitFor, hitRev, itSam, rVar); // Record the other HitRev HSPs if there are any
		else
		{
			hitFor->hit_hsps = hitFor->hit_hsps->next;
			if (hitFor->hit_hsps != NULL)
			{
				hitRev->hit_hsps = rVar->hspTmp;
				rVar->hspTmp = NULL;
				return hitRecord(hitFor, hitRev, itSam, rVar); // Record the other HitFor HSPs if there are any
			}
		}
	}
	
	
	if (hitFor->next != NULL)
	{
		itSam->countHit++;
		itSam->samHits = (SamHitPtr*) saferealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
		itSam->samHits[itSam->countHit-1]->countRec = 0;
		if (hitRev == NULL)
		{
			if (rVar->hitTmp == NULL)
				return hitRecord(hitFor->next, NULL, itSam, rVar);
			else
				return hitRecord(hitFor->next, rVar->hitTmp, itSam, rVar);
		}
		else
			return hitRecord(hitFor->next, hitRev->next, itSam, rVar);
	}
	
	if (rVar->hitTmp != NULL)
	{
		if (!rVar->end)
		{
			hitRev = rVar->hitTmp;
			rVar->end = 1;
		}
	
		if (hitRev->next != NULL)
		{
			for (rVar->tmpHitNb = 0; rVar->tmpHitNb < countHit; rVar->tmpHitNb++)
			{
				if (!strcmp(hitRev->next->hit_def, itSam->samHits[rVar->tmpHitNb]->rsSam[0]->samOut[1]->rname))
					break;
			}
			
			if (rVar->tmpHitNb == countHit)
			{
				itSam->countHit++;
				itSam->samHits = (SamHitPtr*) saferealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
				itSam->samHits[itSam->countHit-1]->countRec = 0;
			}
			return hitRecord(hitFor, hitRev->next, itSam, rVar);
		}
	}
	
	return itSam;	


// for the flag, to substract you can use flag &= ~0x100

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



/*
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
*/


