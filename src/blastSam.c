
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "blastSam.h"
#include "utils.h"
#include "shortRead.h"


typedef struct Cigar
{
	int* str;
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
static int parseDict(char* filename)
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

		fprintf(stdout, "%s", str); // Print the header line in the SAM file
		free(str);
		if (c != EOF)
			fprintf(stdout, "\n");

	} while (c != EOF);

	fclose(reader);

	return 0;
}

#define CSTRMACRO(TAG, FUNCT) { \
		COUNTCHARMACRO(FUNCT); \
		cigar->str[cigar->size-2] = count; \
		cigar->str[cigar->size-1] = TAG;}

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
	int pos = 0;
	int count = 0;
	size_t queryLength = samOut->query->read_len;
	HspPtr hsp = samOut->hsp;
	CigarPtr cigar = (CigarPtr) safeCalloc(1, sizeof(Cigar));

	if (hsp->hsp_query_from > 1) // Soft clipping at the beginning
	{
		cigar->size += 2;
		cigar->str = (int*) safeCalloc(cigar->size, sizeof(int));
		cigar->str[0] = hsp->hsp_query_from - 1;
		cigar->str[1] = 'S';
	}

	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{
		cigar->size += 2;
		cigar->str = (int*) safeRealloc(cigar->str, cigar->size * sizeof(int));

		if (hsp->hsp_hseq[pos] == '-')
			CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-')) // Count the number of insertions

		else if (hsp->hsp_qseq[pos] == '-')
			CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-')) // Count the number of deletions

		else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
			CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])) // Count the number of matches

		else
			CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos])) // Count the number of mismatches

		if (cigar->str[cigar->size-2] >= 100 && cigar->str[cigar->size-1] == 'D')
			cigar->str[cigar->size-1] = 'N'; // If there is more than a hundred deletion at a time, it is considered a skipped region
			
		if (cigar->str[cigar->size-1] != '=')
			cigar->nbDiff += cigar->str[cigar->size-2];
	}

	if ((queryLength - hsp->hsp_query_to) > 0) // Soft clipping at the end
	{
		cigar->size += 2;
		cigar->str = (int*) safeRealloc(cigar->str, cigar->size * sizeof(int));
		cigar->str[cigar->size-2] = queryLength - hsp->hsp_query_to;
		cigar->str[cigar->size-1] = 'S';
	}

	return cigar;
}


/* Record all the alignment hits for a sequence (single end) or a pair of sequences (paired end) */
static IterationSamPtr hitRecord(HitPtr hitFor, HitPtr hitRev, IterationSamPtr itSam, RVarPtr rVar)
{
	int i = 0;
	size_t countHit;
	size_t countRec;
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

	itSam->samHits[countHit-1]->rsSam = (RecordSamPtr*) safeRealloc(itSam->samHits[countHit-1]->rsSam, countRec * sizeof(RecordSamPtr)); // Create and/or append a given Hit record table
	rSam = itSam->samHits[countHit-1]->rsSam[countRec-1] = (RecordSamPtr) safeCalloc(1, sizeof(RecordSam)); // Create a new record

	rSam->samOut[0] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput)); // Create a structure to capture all the info concerning the forward strand
	rSam->samOut[0]->query = rVar->reads[0]; // Record the forward sequence infos

	if (hitRev != NULL)
	{
		rSam->samOut[1] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput));
		rSam->samOut[1]->query = rVar->reads[1]; // Record the reverse sequence infos
	}

	if (hitFor->hit_def == NULL)
	{
		rVar->end = 1;
		rVar->tmpHitNb = countHit;
	}

	if (hitRev == NULL || hitRev->hit_def == NULL) // Single end or the forward strand is mapped on a reference where the reverse strand is not
	{
		if (hitFor->hit_hsps == NULL) // Not mapped
			rSam->score = 0;

		else
		{
			rSam->samOut[0]->rname = shortName(hitFor->hit_def);
			rSam->samOut[0]->hsp = hitFor->hit_hsps;
			rSam->samOut[0]->cigar = cigarStrBuilding(rSam->samOut[0]);
			if (hitRev == NULL)
				rSam->score = 60;
			hitFor->hit_hsps = hitFor->hit_hsps->next;
			if (hitFor->hit_hsps != NULL)
				return hitRecord(hitFor, hitRev, itSam, rVar); // Record the other HSPs if there are any
		}
	}

	else if (rVar->end) // The reverse strand is mapped on a reference where the forward strand is not
	{
		if (rVar->tmpHitNb == countHit)
		{
			rSam->samOut[1]->rname = shortName(hitRev->hit_def);
			rSam->samOut[1]->hsp = hitRev->hit_hsps;
			rSam->samOut[1]->cigar = cigarStrBuilding(rSam->samOut[1]);
			hitRev->hit_hsps = hitRev->hit_hsps->next;
			if (hitRev->hit_hsps != NULL)
				return hitRecord(hitFor, hitRev, itSam, rVar); // Record the other HSPs if there are any
		}
	}

	else if (strcmp(hitFor->hit_def, hitRev->hit_def))
	{
		if (hitRev->next == NULL || !strcmp(hitFor->hit_def, hitRev->next->hit_def))
		{
			itSam->countHit++;
			itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
			itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
		}

		return hitRecord(hitFor, hitRev->next, itSam, rVar);
	}

	else if (strcmp(hitFor->hit_def, hitRev->hit_def) == 0) // Forward and reverse strand are mapped on the same reference
	{
		if (rVar->hspTmp == NULL)
			rVar->hspTmp = hitRev->hit_hsps;

		for (i = 0; i < 2; i++)
		{
			rSam->samOut[i]->hsp = (!i ? hitFor->hit_hsps : hitRev->hit_hsps);
			rSam->samOut[i]->rname = shortName((!i ? hitFor->hit_def : hitRev->hit_def));
			rSam->samOut[i]->cigar = cigarStrBuilding(rSam->samOut[i]);
		}
		rSam->score = 60;
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
		itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
		itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
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
				itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
				itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
			}
			return hitRecord(hitFor, hitRev->next, itSam, rVar);
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

// NOTE: for the flag, to substract you can use flag &= ~0x100

#define SAM_PAIRED 0x1 // Paired end
#define SAM_PROPER_PAIR 0x2 // Read mapped in a proper pair
#define SAM_UNMAP 0x4 // Read unmapped
#define SAM_MUNMAP 0x8 // Mate unmapped
#define SAM_REVERSE 0x10 // Read mapped to the reverse strand
#define SAM_MREVERSE 0x20 // Mate mapped to the reverse strand
#define SAM_READF 0x40 // Read is first in pair
#define SAM_READR 0x80 // Read is last in pair
#define SAM_SECONDARY 0x100 // Not primary alignment
#define SAM_QCFAIL 0x200 // Failed control quality
#define SAM_DUP 0x400 // Optical or PCR duplicate
#define SAM_SUPPLEMENTARY 0x800 // Supplementary alignment

static void printSam(IterationSamPtr itSam)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int invk = 0;
	int pos[2];
	int tlen = 0;
	int p0 = 0;
	int p1 = 0;
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
				rSam = itSam->samHits[i]->rsSam[j];
				if (rSam->samOut[k] == NULL) continue; // Used if single end

				if (rSam->samOut[invk] != NULL) // Paired end
				{
					flag |= SAM_PAIRED | (!k ? SAM_READF : SAM_READR);
					if (rSam->samOut[invk]->rname != NULL)
					{
						flag |= (rSam->samOut[invk]->hsp->hsp_hit_to - rSam->samOut[invk]->hsp->hsp_hit_from < 0 ? SAM_MREVERSE : 0);
						pos[invk] = (flag & SAM_MREVERSE ? rSam->samOut[invk]->hsp->hsp_hit_to : rSam->samOut[invk]->hsp->hsp_hit_from);
						if (rSam->samOut[k]->rname != NULL && !strcmp(rSam->samOut[k]->rname, rSam->samOut[invk]->rname))
						{
							rnext = "=";
							flag |= SAM_PROPER_PAIR;
							flag |= (rSam->samOut[k]->hsp->hsp_hit_to - rSam->samOut[k]->hsp->hsp_hit_from < 0 ? SAM_REVERSE : 0);
							p0 = (flag & SAM_REVERSE ? rSam->samOut[k]->hsp->hsp_hit_to + rSam->samOut[k]->hsp->hsp_align_len : rSam->samOut[k]->hsp->hsp_hit_from);
							p1 = (flag & SAM_MREVERSE ? rSam->samOut[invk]->hsp->hsp_hit_to + rSam->samOut[invk]->hsp->hsp_align_len : rSam->samOut[invk]->hsp->hsp_hit_from);
							tlen = p1 - p0;
						}
						else
							rnext = rSam->samOut[invk]->rname;
						
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
					fprintf(stdout, "%s\t%d\t*\t0\t%d\t*\t%s\t%d\t0\t%s\t%s\n", rSam->samOut[k]->query->name, flag, rSam->score, rnext, pos[invk], rSam->samOut[k]->query->seq, rSam->samOut[k]->query->qual);
				}
				else
				{
					flag |= (rSam->samOut[k]->hsp->hsp_hit_to - rSam->samOut[k]->hsp->hsp_hit_from < 0 ? SAM_REVERSE : 0);
					pos[k] = (flag & SAM_REVERSE ? rSam->samOut[k]->hsp->hsp_hit_to : rSam->samOut[k]->hsp->hsp_hit_from);
					fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t", rSam->samOut[k]->query->name, flag, rSam->samOut[k]->rname, pos[k], rSam->score);
					if (flag & SAM_REVERSE)
					{
						for (l = rSam->samOut[k]->cigar->size - 2; l >= 0; l -= 2)
							fprintf(stdout,"%d%c", rSam->samOut[k]->cigar->str[l], rSam->samOut[k]->cigar->str[l+1]);
						seq = revStr(rSam->samOut[k]->query->seq);
					}
					else
					{
						for (l = 0; l < (rSam->samOut[k]->cigar->size - 1); l += 2)
							fprintf(stdout,"%d%c", rSam->samOut[k]->cigar->str[l], rSam->samOut[k]->cigar->str[l+1]);
						seq = rSam->samOut[k]->query->seq;
					}
					fprintf(stdout, "\t%s\t%d\t%d\t%s\t", rnext, pos[invk], tlen, seq);
					if (flag & SAM_REVERSE)
					{
						free(seq);
						for (l = (strlen(rSam->samOut[k]->query->qual) - 1); l >= 0; l--)
							fprintf(stdout, "%c", rSam->samOut[k]->query->qual[l]);
					}
					else
						fprintf(stdout, "%s", rSam->samOut[k]->query->qual);
						
					fprintf(stdout, "\tNM:i:%d\n", rSam->samOut[k]->cigar->nbDiff);
				}
			}
		}
	}
}

static IterationSamPtr iterationRecord(xmlTextReaderPtr reader, gzFile fp, gzFile fp2, int inter)
{
	IterationSamPtr itSam = NULL;
	IterationPtr itFor = NULL;
	IterationPtr itRev = NULL;
	RVarPtr rVar = (RVarPtr) safeCalloc(1, sizeof(RVar));

	if (fp2 != NULL || inter) // Paired end
	{
		itFor = parseIteration(reader);
		itRev = parseIteration(reader);
		rVar->hitTmp = itRev->iteration_hits;

		rVar->reads[0] = shortReadNext(fp);
		if (inter)
			rVar->reads[1] = shortReadNext(fp);
		else
			rVar->reads[1] = shortReadNext(fp2);

		itSam = hitRecord(itFor->iteration_hits, itRev->iteration_hits, itSam, rVar);
		printSam(itSam);

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
		printSam(itSam);
		deallocIteration(itFor);
		shortReadFree(rVar->reads[0]);
	}

	free(rVar);

	return itSam;
}


static void deallocItSam(IterationSamPtr itSam)
{
	int i = 0;
	int j = 0;
	int k = 0;

	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				if (itSam->samHits[i]->rsSam[j]->samOut[k] != NULL)
				{
					if (itSam->samHits[i]->rsSam[j]->samOut[k]->cigar != NULL)
					{
						free(itSam->samHits[i]->rsSam[j]->samOut[k]->cigar->str);
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

int blastToSam(int argc, char** argv)
{
	xmlTextReaderPtr reader;
	gzFile fp = NULL;
	gzFile fp2 = NULL;
	BlastOutputPtr blastOP = NULL;
	IterationSamPtr itSam = NULL;
	int inter = 0; // TODO: Put it in option -> get_longopt/getopt (look at Jennifer's code)

	reader = safeXmlNewTextReaderFilename(argv[1]);

	safeXmlTextReaderRead(reader);

	if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
		ERROR("The document is not a Blast output\n", 1)

	safeXmlTextReaderRead(reader);

	blastOP = parseBlastOutput(reader);

	if (parseDict(argv[2]) == 1)
		ERROR("Error while reading the index base", 1)

	fp = safeGzOpen(argv[3], "r");

	if (argc == 5 && !inter)
		fp2 = safeGzOpen(argv[4], "r");

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		itSam = iterationRecord(reader, fp, fp2, inter);
		if (itSam != NULL)
			deallocItSam(itSam);
	}

	gzclose(fp);
	if (fp2 != NULL)
		gzclose(fp2);

	deallocBlastOutput(blastOP);

	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	return 0;
}


