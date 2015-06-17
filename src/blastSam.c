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
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "blastSam.h"
#include "utils.h"
#include "shortRead.h"
#include "debug.h"

/************************************************************************************/
/*	TODO	*/
/************************************************************************************/
typedef struct RecordVariables // Temp structure used in hitRecord()
{
	int end; // Marks the fact that all the first in pair hits have been recorded
	
	// These variables are used to record second in pair read's alignment on references where the first in pair is not mapped
	int tmpHitNb; // Number of references where the second in pair read is mapped
	HitPtr hitTmp; // Temp pointer to the first reference where the second in pair read is mapped 
	
	// Used to enable the record of all first in pair alignments with all second in pair alignments
	HspPtr hspTmp; // Temp pointer to the first HSP of the second in pair read
	
	ShortReadPtr reads[2]; // reads infos from the fastQ. 0: first in pair; 1: second in pair
} RVar, *RVarPtr;

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
		if (hitFirst->hit_hsps != NULL) // Mapped
		{
			rSam->samOut[0]->rname = safeStrdup(hitFirst->hit_def);
			rSam->samOut[0]->hsp = hitFirst->hit_hsps;
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
		}
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
				if (!strcmp(hitSec->next->hit_def, itSam->samHits[rVar->tmpHitNb]->rsSam[0]->samOut[1]->rname)) break;

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


/************************************************************************************/
/*	TODO	*/
/************************************************************************************/
static IterationSamPtr iterationRecord(xmlTextReaderPtr reader, gzFile fp, gzFile fp2, AppParamPtr app)
{
	IterationSamPtr itSam = NULL;
	IterationPtr itFirst = NULL;
	IterationPtr itSec = NULL;
	RVarPtr rVar = (RVarPtr) safeCalloc(1, sizeof(RVar));

	if (fp2 != NULL || app->inter) // Paired end
	{
		itFirst = parseIteration(reader);
		itSec = parseIteration(reader);
		rVar->hitTmp = itSec->iteration_hits;

		rVar->reads[0] = shortReadNext(fp);
		if (app->inter)
			rVar->reads[1] = shortReadNext(fp);
		else
			rVar->reads[1] = shortReadNext(fp2);

		itSam = hitRecord(itFirst->iteration_hits, itSec->iteration_hits, itSam, rVar);
		printSam(itSam, app);

		deallocIteration(itFirst);
		deallocIteration(itSec);
		shortReadFree(rVar->reads[0]);
		shortReadFree(rVar->reads[1]);
	}

	else // Single end
	{
		itFirst = parseIteration(reader);
		rVar->reads[0] = shortReadNext(fp);
		itSam = hitRecord(itFirst->iteration_hits, NULL, itSam, rVar);
		printSam(itSam, app);
		deallocIteration(itFirst);
		shortReadFree(rVar->reads[0]);
	}

	free(rVar);
	return itSam;
}


/************************************************************************************/
/*	TODO	*/
/************************************************************************************/
static void deallocItSam(IterationSamPtr itSam)
{
	int i = 0, j = 0, k = 0;

	for (i = 0; i < itSam->countHit; i++)
	{
		for (j = 0; j < itSam->samHits[i]->countRec; j++)
		{
			for (k = 0; k < 2; k++)
			{
				if (itSam->samHits[i]->rsSam[j]->samOut[k] != NULL)
				{
					if (itSam->samHits[i]->rsSam[j]->samOut[k]->rname != NULL)
						free(itSam->samHits[i]->rsSam[j]->samOut[k]->rname);
					free(itSam->samHits[i]->rsSam[j]->samOut[k]);
				}
			}
			free(itSam->samHits[i]->rsSam[j]);
		}
		free(itSam->samHits[i]);
	}
	free(itSam);
}


/************************************************************************************/
/*	TODO	*/
/************************************************************************************/
static char* readGroupID(char* readGroup)
{
	char* rgID = NULL;
	char* str = strstr(readGroup, "\tID:");
	str += 4;
	size_t i, str_length = strlen(str);
	int c = 0;
	
	for (i = 0; i <= str_length && c != '\t'; i++)
		c = str[i];

	rgID = safeMalloc(i);
	memcpy(rgID, str, i-1);
	rgID[i-1] = '\0';
	return rgID;
}

int blastToSam(AppParamPtr app)
{
	xmlTextReaderPtr reader;
	gzFile fp = NULL;
	gzFile fp2 = NULL;
	BlastOutputPtr blastOP = NULL;
	IterationSamPtr itSam = NULL;

	reader = safeXmlNewTextReaderFilename(app->blastOut);

	safeXmlTextReaderRead(reader);

	if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
		ERROR("The document is not a Blast output\n", 1);

	safeXmlTextReaderRead(reader);

	blastOP = parseBlastOutput(reader);

	if (samHead(app) == 1)
		ERROR("Error while printing the Sam header\n", 1);

	fp = safeGzOpen(app->fastq1, "r");

	if (app->fastq2 != NULL)
		fp2 = safeGzOpen(app->fastq2, "r");

	if (app->readGroup != NULL)
		app->readGroupID = readGroupID(app->readGroup);

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		itSam = iterationRecord(reader, fp, fp2, app);
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


