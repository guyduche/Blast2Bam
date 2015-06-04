
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "blastSam.h"
#include "utils.h"
#include "shortRead.h"

typedef struct CigarElement
{
	int count;
	char symbol;
} CigarElement,*CigarElementPtr;

typedef struct Cigar
{
	CigarElementPtr* elements;
	int nbDiff;
	size_t size; // size of the CIGAR string
} Cigar, *CigarPtr;

typedef struct SamOutput
{
	ShortReadPtr query; // Sequence infos
	HspPtr hsp; // HSP infos
	CigarPtr cigar; // CIGAR string
	char* rname; // Ref name
} SamOutput, *SamOutputPtr;

typedef struct RecordSam
{
	SamOutputPtr samOut[2]; // 0: first read; 1: mate
	unsigned int score; // Mapping quality
} RecordSam, *RecordSamPtr;

typedef struct SamHit
{
	RecordSamPtr* rsSam; // List of all the paired records of a reference hit
	size_t countRec; // Number of paired records
	size_t countHSPsec; // Number of HSP found for the second read (useful for option -W)
} SamHit, *SamHitPtr;

typedef struct IterationSam
{
	SamHitPtr* samHits; // List of all the reference hits
	size_t countHit; // Number of hits (number of references on which an alignment has been found for the sequence)
} IterationSam, *IterationSamPtr;

typedef struct RecordVariables // Temp structure used in hitRecord()
{
	int end;
	int tmpHitNb;
	HitPtr hitTmp;
	HspPtr hspTmp;
	ShortReadPtr reads[2]; // sequences infos
} RVar, *RVarPtr;


/* Print SAM header */
static void sq_line(char* filename)
{
	FILE* reader;
	char* str = NULL;
	int c;
	int countSpace = 0;
	size_t lenStr = 0;

	reader = safeFOpen(filename, "r");

	do
	{
		if (countSpace > 2 || !countSpace) // Skip the first line of the file and the end of the other lines
			do
			{
				c = fgetc(reader);
			} while (c != '\n' && c != EOF);

		lenStr = 0;
		countSpace = 0;
		str = NULL;
		c = fgetc(reader);
		lenStr++;

		while (c != '\n' && c != EOF && countSpace <= 2) // Keep only the reference name and its length
		{
			if (c == '\t')
				countSpace++;

			str = (char*) safeRealloc(str, lenStr+1);
			str[lenStr-1] = (char) c;
			c = fgetc(reader);
			lenStr++;
		}

		str = (char*) safeRealloc(str, lenStr);
		str[lenStr-1] = '\0';

		fprintf(stdout, "%s", str); // Print the SQ line in the SAM file
		free(str);
		if (c != EOF)
			fprintf(stdout, "\n");

	} while (c != EOF);

	fclose(reader);
}

static int rg_line(char* readGroup)
{	
	if (strstr(readGroup, "@RG") != readGroup) return 1;
	if (strstr(readGroup, "\tID:") == NULL) return 1;
	fprintf(stdout, "%s\n", readGroup);
	return 0;
}

static int samHead(AppParamPtr app)
{
	sq_line(app->db);
	if (app->readGroup != NULL)
		if (rg_line(app->readGroup) == 1) return 1;
	return 0;
}

static char* readGroupID(char* readGroup)
{
	size_t i;
	int c = 0;
	char* rgID = NULL;
	char* str = strstr(readGroup, "\tID:");
	str += 4;
	
	for (i = 0; i <= strlen(str) && c != '\t'; i++)
		c = str[i];

	rgID = safeMalloc(i);
	memcpy(rgID, str, i-1);
	rgID[i-1] = '\0';
	return rgID;
}

#define CSTRMACRO(TAG, FUNCT) { \
		COUNTCHARMACRO(FUNCT); \
		cigar->elements[cigar->size-1]->count = count; \
		cigar->elements[cigar->size-1]->symbol = TAG;}

#define COUNTCHARMACRO(FUNCT) { \
		count = 1; \
		pos++; \
		while (pos < (hsp->hsp_align_len) && FUNCT) \
		{ \
			count++; \
			pos++; \
		} \
		pos--;}

/* Build the CIGAR string of the query */
static CigarPtr cigarStrBuilding(SamOutputPtr samOut)
{
	int pos = 0, count = 0;
	size_t queryLength = samOut->query->read_len;
	HspPtr hsp = samOut->hsp;
	CigarPtr cigar = (CigarPtr) safeCalloc(1, sizeof(Cigar));

	if (hsp->hsp_query_from > 1) // Soft clipping at the beginning
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeCalloc(cigar->size, sizeof(CigarElementPtr));
		cigar->elements[0] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));
		cigar->elements[0]->count = hsp->hsp_query_from - 1;
		cigar->elements[0]->symbol = 'S';
	}

	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr));
		cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));

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
			cigar->nbDiff += cigar->elements[cigar->size-1]->count;
	}

	if ((queryLength - hsp->hsp_query_to) > 0) // Soft clipping at the end
	{
		cigar->size++;
		cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr));
		cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));
		cigar->elements[cigar->size-1]->count = queryLength - hsp->hsp_query_to;
		cigar->elements[cigar->size-1]->symbol = 'S';
	}

	return cigar;
}

/* Record all the alignment hits for a sequence (single end) or a pair of sequences (paired end) */
static IterationSamPtr hitRecord(HitPtr hitFirst, HitPtr hitSec, IterationSamPtr itSam, RVarPtr rVar)
{
	int i = 0;
	size_t countHit, countRec;
	RecordSamPtr rSam;

	if (itSam == NULL) // Create the main structure on the first call
	{
		itSam = (IterationSamPtr) safeCalloc(1, sizeof(IterationSam));
		itSam->countHit++;
		itSam->samHits = (SamHitPtr*) safeMalloc(sizeof(SamHitPtr));
		itSam->samHits[0] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
	}

	countHit = itSam->countHit;
	countRec = ++(itSam->samHits[countHit-1]->countRec);

	itSam->samHits[countHit-1]->rsSam = (RecordSamPtr*) safeRealloc(itSam->samHits[countHit-1]->rsSam, countRec * sizeof(RecordSamPtr)); // Create or append a given Hit record table
	rSam = itSam->samHits[countHit-1]->rsSam[countRec-1] = (RecordSamPtr) safeCalloc(1, sizeof(RecordSam)); // Create a new record

	rSam->samOut[0] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput)); // Create a structure to capture all the info concerning the first read
	rSam->samOut[0]->query = rVar->reads[0]; // Record the first read infos

	if (hitSec != NULL)
	{
		rSam->samOut[1] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput));
		rSam->samOut[1]->query = rVar->reads[1]; // Record the second read infos
	}

	if (hitFirst->hit_def == NULL)
	{
		rVar->end = 1;
		rVar->tmpHitNb = countHit;
	}

	if (hitSec == NULL || hitSec->hit_def == NULL) // Single end or the first read is mapped on a reference where his mate is not
	{
		if (hitFirst->hit_hsps == NULL) // Not mapped
			rSam->score = 0;

		else
		{
			rSam->samOut[0]->rname = safeStrdup(hitFirst->hit_def);
			rSam->samOut[0]->hsp = hitFirst->hit_hsps;
			rSam->samOut[0]->cigar = cigarStrBuilding(rSam->samOut[0]);
			if (hitSec == NULL)
				rSam->score = 60;
			hitFirst->hit_hsps = hitFirst->hit_hsps->next;
			if (hitFirst->hit_hsps != NULL)
				return hitRecord(hitFirst, hitSec, itSam, rVar); // Record the other HSPs if there are any
		}
	}

	else if (rVar->end) // The second read is mapped on a reference where the first read is not
	{
		if (rVar->tmpHitNb == countHit)
		{
			itSam->samHits[countHit-1]->countHSPsec++;
			rSam->samOut[1]->rname = safeStrdup(hitSec->hit_def);
			rSam->samOut[1]->hsp = hitSec->hit_hsps;
			rSam->samOut[1]->cigar = cigarStrBuilding(rSam->samOut[1]);
			hitSec->hit_hsps = hitSec->hit_hsps->next;
			if (hitSec->hit_hsps != NULL)
				return hitRecord(hitFirst, hitSec, itSam, rVar); // Record the other HSPs if there are any
		}
	}

	else if (strcmp(hitFirst->hit_def, hitSec->hit_def))
	{
		if (hitSec->next == NULL || !strcmp(hitFirst->hit_def, hitSec->next->hit_def))
		{
			itSam->countHit++;
			itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
			itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
		}

		return hitRecord(hitFirst, hitSec->next, itSam, rVar);
	}

	else if (strcmp(hitFirst->hit_def, hitSec->hit_def) == 0) // Both reads are mapped on the same reference
	{
		if (rVar->hspTmp == NULL)
			rVar->hspTmp = hitSec->hit_hsps;

		itSam->samHits[countHit-1]->countHSPsec++;
		for (i = 0; i < 2; i++)
		{
			rSam->samOut[i]->hsp = (!i ? hitFirst->hit_hsps : hitSec->hit_hsps);
			rSam->samOut[i]->rname = safeStrdup((!i ? hitFirst->hit_def : hitSec->hit_def));
			rSam->samOut[i]->cigar = cigarStrBuilding(rSam->samOut[i]);
		}
		rSam->score = 60;
		hitSec->hit_hsps = hitSec->hit_hsps->next;
		if (hitSec->hit_hsps != NULL)
			return hitRecord(hitFirst, hitSec, itSam, rVar); // Record the other HitSec HSPs if there are any
		else
		{
			hitFirst->hit_hsps = hitFirst->hit_hsps->next;
			if (hitFirst->hit_hsps != NULL)
			{
				hitSec->hit_hsps = rVar->hspTmp;
				rVar->hspTmp = NULL;
				itSam->samHits[countHit-1]->countHSPsec = 0;
				return hitRecord(hitFirst, hitSec, itSam, rVar); // Record the other HitFirst HSPs if there are any
			}
		}
	}
	
	if (hitFirst->next != NULL)
	{
		itSam->countHit++;
		itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
		itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
		if (hitSec == NULL)
		{
			if (rVar->hitTmp == NULL)
				return hitRecord(hitFirst->next, NULL, itSam, rVar);
			else
				return hitRecord(hitFirst->next, rVar->hitTmp, itSam, rVar);
		}
		else
			return hitRecord(hitFirst->next, hitSec->next, itSam, rVar);
	}

	if (rVar->hitTmp != NULL)
	{
		if (!rVar->end)
		{
			hitSec = rVar->hitTmp;
			rVar->end = 1;
		}

		if (hitSec->next != NULL)
		{
			for (rVar->tmpHitNb = 0; rVar->tmpHitNb < countHit; rVar->tmpHitNb++)
			{
				if (!strcmp(hitSec->next->hit_def, itSam->samHits[rVar->tmpHitNb]->rsSam[0]->samOut[1]->rname))
					break;
			}

			if (rVar->tmpHitNb == countHit)
			{
				itSam->countHit++;
				itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
				itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
			}
			return hitRecord(hitFirst, hitSec->next, itSam, rVar);
		}
	}
	return itSam;
}

static char* revStr(char* oldStr)
{
	size_t i = 0;
	int l = strlen(oldStr);
	char* newStr = safeCalloc(l+1, sizeof(char));

	for (--l; l >= 0; l--, i++)
	{
		switch(oldStr[l])
		{
			case 'A': newStr[i] = 'T'; break;
			case 'T': newStr[i] = 'A'; break;
			case 'C': newStr[i] = 'G'; break;
			case 'G': newStr[i] = 'C'; break;
			default: newStr[i] = oldStr[l]; break;
		}
	}
	return newStr;
}

static int firstPosRef(char* rname)
{
	int i, c = 0;
	
	for (i = 0; i < strlen(rname) && c != ':'; i++)
		c = rname[i];
	
	return strtol(rname + i, NULL, 10);
}

// NOTE: for the flag, to substract you can use flag &= ~0x100

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

static void printSam(IterationSamPtr itSam, AppParamPtr app, char* rgID)
{
	int i = 0, j = 0, k = 0, l = 0;
	int invk = 0, tlen = 0, p0 = 0, p1 = 0, doNotPrint = 0, countUnprint = 0;
	long int pos[2];
	unsigned int flag = 0;
	char* rnext = NULL;
	char* seq = NULL;
	RecordSamPtr rSam;

	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				flag = 0;
				invk = (k ? 0 : 1);
				if (!k) doNotPrint = 0;
				rSam = itSam->samHits[i]->rsSam[j];
				if (rSam->samOut[k] == NULL || doNotPrint) continue;
				if (app->minLen != 0 && !k)
				{
					if (rSam->samOut[0]->hsp != NULL && rSam->samOut[0]->hsp->hsp_align_len > app->minLen)
					{
						if (rSam->samOut[1] != NULL && rSam->samOut[1]->hsp != NULL && rSam->samOut[1]->hsp->hsp_align_len < app->minLen)
						{
							if (countUnprint != (itSam->samHits[i]->countRec - itSam->samHits[i]->countHSPsec))
								{doNotPrint = 1; countUnprint++; continue;}	
							else
								rSam->samOut[1]->hsp = NULL;
						}
					}
					else
					{
						if (rSam->samOut[1] != NULL && rSam->samOut[1]->hsp != NULL && rSam->samOut[1]->hsp->hsp_align_len > app->minLen)
						{
							if (countUnprint != (itSam->samHits[i]->countRec - itSam->samHits[i]->countHSPsec))
								{doNotPrint = 1; countUnprint++; continue;}
							else
								rSam->samOut[0]->hsp = NULL;
						}
						else
						{
							if (countUnprint != (itSam->samHits[i]->countRec -1))
								{doNotPrint = 1; countUnprint++; continue;}
							else
							{
								rSam->samOut[0]->hsp = NULL;
								if (rSam->samOut[1] != NULL)
									rSam->samOut[1]->hsp = NULL;
							}
						}
					}
				}

				if (rSam->samOut[invk] != NULL) // Paired end
				{
					flag |= SAM_PAIRED | (!k ? SAM_READF : SAM_READL);
					if (rSam->samOut[invk]->hsp != NULL)
					{
						flag |= (rSam->samOut[invk]->hsp->hsp_hit_to - rSam->samOut[invk]->hsp->hsp_hit_from < 0 ? SAM_MREVERSE : 0);
						pos[invk] = (flag & SAM_MREVERSE ? rSam->samOut[invk]->hsp->hsp_hit_to : rSam->samOut[invk]->hsp->hsp_hit_from);
						if (app->posOnChr) pos[invk] += firstPosRef(rSam->samOut[invk]->rname);
						if (rSam->samOut[k]->hsp != NULL && !strcmp(rSam->samOut[k]->rname, rSam->samOut[invk]->rname))
						{
							rnext = "=";
							flag |= SAM_PROPER_PAIR;
							flag |= (rSam->samOut[k]->hsp->hsp_hit_to - rSam->samOut[k]->hsp->hsp_hit_from < 0 ? SAM_REVERSE : 0);
							p0 = (flag & SAM_REVERSE ? rSam->samOut[k]->hsp->hsp_hit_to + rSam->samOut[k]->hsp->hsp_align_len : rSam->samOut[k]->hsp->hsp_hit_from);
							p1 = (flag & SAM_MREVERSE ? rSam->samOut[invk]->hsp->hsp_hit_to + rSam->samOut[invk]->hsp->hsp_align_len : rSam->samOut[invk]->hsp->hsp_hit_from);
							tlen = p1 - p0;
						}
						else
							rnext = shortName(rSam->samOut[invk]->rname);
					}
					else
					{
						flag |= SAM_MUNMAP;
						rnext = "*";
						pos[invk] = 0;
					}
				}
				else // Single end
				{
					rnext = "*";
					pos[invk] = 0;
				}

				if (rSam->samOut[k]->hsp == NULL)
				{
					flag |= SAM_UNMAP;
					fprintf(stdout, "%s\t%d\t*\t0\t%d\t*\t%s\t%ld\t0\t%s\t%s", rSam->samOut[k]->query->name, flag, rSam->score, rnext, pos[invk], rSam->samOut[k]->query->seq, rSam->samOut[k]->query->qual);
				}
				else
				{
					flag |= (rSam->samOut[k]->hsp->hsp_hit_to - rSam->samOut[k]->hsp->hsp_hit_from < 0 ? SAM_REVERSE : 0);
					pos[k] = (flag & SAM_REVERSE ? rSam->samOut[k]->hsp->hsp_hit_to : rSam->samOut[k]->hsp->hsp_hit_from);
					if (app->posOnChr) pos[k] += firstPosRef(rSam->samOut[k]->rname);
					fprintf(stdout, "%s\t%d\t%s\t%ld\t%d\t", rSam->samOut[k]->query->name, flag, shortName(rSam->samOut[k]->rname), pos[k], rSam->score);
					if (flag & SAM_REVERSE)
					{
						for (l = rSam->samOut[k]->cigar->size - 1; l >= 0; l--)
							fprintf(stdout,"%d%c", rSam->samOut[k]->cigar->elements[l]->count, rSam->samOut[k]->cigar->elements[l]->symbol);
						seq = revStr(rSam->samOut[k]->query->seq);
					}
					else
					{
						for (l = 0; l < rSam->samOut[k]->cigar->size; l++)
							fprintf(stdout,"%d%c", rSam->samOut[k]->cigar->elements[l]->count, rSam->samOut[k]->cigar->elements[l]->symbol);
						seq = rSam->samOut[k]->query->seq;
					}
					fprintf(stdout, "\t%s\t%ld\t%d\t%s\t", rnext, pos[invk], tlen, seq);
					if (flag & SAM_REVERSE)
					{
						free(seq);
						for (l = (strlen(rSam->samOut[k]->query->qual) - 1); l >= 0; l--)
							fprintf(stdout, "%c", rSam->samOut[k]->query->qual[l]);
					}
					else
						fprintf(stdout, "%s", rSam->samOut[k]->query->qual);

					fprintf(stdout, "\tNM:i:%d", rSam->samOut[k]->cigar->nbDiff);
				}
				if (rgID != NULL)
					fprintf(stdout, "\tRG:Z:%s", rgID);
				fprintf(stdout, "\n");
			}
		}
		countUnprint = 0;
	}
}

static IterationSamPtr iterationRecord(xmlTextReaderPtr reader, gzFile fp, gzFile fp2, AppParamPtr app, char* rgID)
{
	IterationSamPtr itSam = NULL;
	IterationPtr itFor = NULL;
	IterationPtr itRev = NULL;
	RVarPtr rVar = (RVarPtr) safeCalloc(1, sizeof(RVar));

	if (fp2 != NULL || app->inter) // Paired end
	{
		itFor = parseIteration(reader);
		itRev = parseIteration(reader);
		rVar->hitTmp = itRev->iteration_hits;

		rVar->reads[0] = shortReadNext(fp);
		if (app->inter)
			rVar->reads[1] = shortReadNext(fp);
		else
			rVar->reads[1] = shortReadNext(fp2);

		itSam = hitRecord(itFor->iteration_hits, itRev->iteration_hits, itSam, rVar);
		printSam(itSam, app, rgID);

		deallocIteration(itFor);
		deallocIteration(itRev);
		shortReadFree(rVar->reads[0]);
		shortReadFree(rVar->reads[1]);
	}

	else // Single end
	{
		itFor = parseIteration(reader);
		rVar->reads[0] = shortReadNext(fp);
		itSam = hitRecord(itFor->iteration_hits, NULL, itSam, rVar);
		printSam(itSam, app, rgID);
		deallocIteration(itFor);
		shortReadFree(rVar->reads[0]);
	}

	free(rVar);
	return itSam;
}


static void deallocItSam(IterationSamPtr itSam)
{
	int i = 0, j = 0, k = 0, l = 0;

	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				if (itSam->samHits[i]->rsSam[j]->samOut[k] != NULL)
				{
					if (itSam->samHits[i]->rsSam[j]->samOut[k]->rname != NULL)
					{
						free(itSam->samHits[i]->rsSam[j]->samOut[k]->rname);
						for (l = 0; l < itSam->samHits[i]->rsSam[j]->samOut[k]->cigar->size; l++)
							free(itSam->samHits[i]->rsSam[j]->samOut[k]->cigar->elements[l]);
						free(itSam->samHits[i]->rsSam[j]->samOut[k]->cigar->elements);
						free(itSam->samHits[i]->rsSam[j]->samOut[k]->cigar);
					}
					free(itSam->samHits[i]->rsSam[j]->samOut[k]);
				}
			}
			free(itSam->samHits[i]->rsSam[j]);
		}
		free(itSam->samHits[i]);
	}
	free(itSam);
}

int blastToSam(AppParamPtr app)
{
	xmlTextReaderPtr reader;
	gzFile fp = NULL;
	gzFile fp2 = NULL;
	BlastOutputPtr blastOP = NULL;
	IterationSamPtr itSam = NULL;
	char* rgID = NULL;

	reader = safeXmlNewTextReaderFilename(app->blastOut);

	safeXmlTextReaderRead(reader);

	if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
		ERROR("The document is not a Blast output\n", 1)

	safeXmlTextReaderRead(reader);

	blastOP = parseBlastOutput(reader);

	if (samHead(app) == 1)
		ERROR("Error while printing the Sam header\n", 1)

	fp = safeGzOpen(app->fastq1, "r");

	if (app->fastq2 != NULL)
		fp2 = safeGzOpen(app->fastq2, "r");

	if (app->readGroup != NULL)
		rgID = readGroupID(app->readGroup);

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		itSam = iterationRecord(reader, fp, fp2, app, rgID);
		if (itSam != NULL)
			deallocItSam(itSam);
	}

	gzclose(fp);
	if (fp2 != NULL)
		gzclose(fp2);

	deallocBlastOutput(blastOP);
	if (rgID != NULL) free(rgID);
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	return 0;
}


