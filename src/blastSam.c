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
#include <string.h>
#include "blastSam.h"
#include "utils.h"

/************************************************************************************/
/*  Record the alignment hits of a read or a pair of reads (single or paired end)   */
/************************************************************************************/
// Temp structure used in hitRecord()
typedef struct RecordVariables
{
    int end;                    // Marks the fact that all the first in pair hits have been recorded

                                // These variables are used to record second in pair read's alignment on references where the first in pair is not mapped
    int tmpHitNb;               // Number of references where the second in pair read is mapped
    HitPtr hitTmp;              // Temp pointer to the first reference where the second in pair read is mapped

    HspPtr hspTmp;              // Temp pointer to the first HSP of the second in pair read
                                // Used to enable the record of all first in pair alignments with all second in pair alignments

    ShortReadPtr reads[2];      // Reads infos from the fastQ. 0: first in pair; 1: second in pair
} RVar, *RVarPtr;


static IterationSamPtr hitRecord(
    HitPtr hitFirst,            // List of Blast hits of the first in pair
    HitPtr hitSec,              // List of Blast hits of the second in pair (NULL in case of single end)
    IterationSamPtr itSam,      // Output
    RVarPtr rVar)               // HitRecord() is recursive, we keep the state of the function in this variable
{
    int i = 0;
    size_t countHit, countRec;
    SamOutputPtr samOut[2] = {NULL, NULL};

    // Create the main structure on the first call
    if (itSam == NULL)
    {
        itSam = (IterationSamPtr) safeCalloc(1, sizeof(IterationSam));      // Main structure
        itSam->countHit++;                                                  // Even if the read is not mapped, it still needs to be recorded
        itSam->samHits = (SamHitPtr*) safeMalloc(sizeof(SamHitPtr));        // Array containing the all the reference hits of the read(s)
        itSam->samHits[0] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));      // Structure containing all the records of the read(s) on a single reference
    }

    countHit = itSam->countHit;
    countRec = ++(itSam->samHits[countHit-1]->countRec);
    itSam->countRecGlobal++;

    itSam->samHits[countHit-1]->rsSam = (RecordSamPtr*) safeRealloc(itSam->samHits[countHit-1]->rsSam, countRec * sizeof(RecordSamPtr));    // Create or append a given Hit record array
    itSam->samHits[countHit-1]->rsSam[countRec-1] = (RecordSamPtr) safeCalloc(1, sizeof(RecordSam));                                        // Create a new record

    samOut[0] = itSam->samHits[countHit-1]->rsSam[countRec-1]->samOut[0] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput));     // Create a structure to capture all the info concerning the first read
    samOut[0]->query = rVar->reads[0];                                                                                          // Record the first read infos

    // Only if paired end
    if (hitSec != NULL)
    {
        samOut[1] = itSam->samHits[countHit-1]->rsSam[countRec-1]->samOut[1] = (SamOutputPtr) safeCalloc(1, sizeof(SamOutput)); // Create a structure to capture all the info concerning the second read
        samOut[1]->query = rVar->reads[1];                                                                                      // Record the second read infos
    }

    // Go there if the first read is unmapped
    if (hitFirst->hit_def == NULL)
    {
        rVar->end = 1;      // Used to enable the possibility that the second read can be mapped on a reference where the first is not
        rVar->tmpHitNb = countHit;
    }

    // Single end or the second read is unmapped
    if (hitSec == NULL || hitSec->hit_def == NULL)
    {
        // The first read is mapped
        if (hitFirst->hit_hsps != NULL)
        {
            samOut[0]->rname = safeStrdup(hitFirst->hit_def);       // Get the reference name
            samOut[0]->hsp = hitFirst->hit_hsps;                    // Get the alignment infos
            hitFirst->hit_hsps = hitFirst->hit_hsps->next;          // Go to the next alignment on the same reference
            if (hitFirst->hit_hsps != NULL)
                return hitRecord(hitFirst, hitSec, itSam, rVar);    // Record the other HSPs on this reference if there are any
        }
    }

    // The second read is mapped on a reference where the first read is not
    else if (rVar->end)
    {
        // Go there only if the second read HSPs of the current reference have not been recorded
        if (rVar->tmpHitNb == countHit)
        {
            itSam->samHits[countHit-1]->countHSPsec++;              // Count the number of HSPs of the second read on the current reference. Useful for option -W
            samOut[1]->rname = safeStrdup(hitSec->hit_def);         // Get the reference name
            samOut[1]->hsp = hitSec->hit_hsps;                      // Get the alignment infos
            hitSec->hit_hsps = hitSec->hit_hsps->next;              // Go to the next alignment on the same reference
            if (hitSec->hit_hsps != NULL)
                return hitRecord(hitFirst, hitSec, itSam, rVar);    // Record the other HSPs on this reference if there are any
        }
    }

    // Both reads are mapped but not on the same reference
    else if (strcmp(hitFirst->hit_def, hitSec->hit_def))
    {
        // Go through all the reference hits of the second read in order to find a match with the current reference of the first read
        // If the next reference hit of the second read is a match, a new Hit structure is created and, at the next function call, it will go through the next section
        // If there is no next reference hit, a new Hit structure is created and, at the next function call, the second read will be recorded as unmapped
        if (hitSec->next == NULL || !strcmp(hitFirst->hit_def, hitSec->next->hit_def))
        {
            itSam->countHit++;
            itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
            itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
        }

        return hitRecord(hitFirst, hitSec->next, itSam, rVar);
    }

    // Both reads are mapped on the same reference
    else if (strcmp(hitFirst->hit_def, hitSec->hit_def) == 0)
    {
        if (rVar->hspTmp == NULL)
            rVar->hspTmp = hitSec->hit_hsps;                        // Keep a pointer on the first HSP of the second read

        itSam->samHits[countHit-1]->countHSPsec++;                  // Number of HSPs of the second read (-W)
        for (i = 0; i < 2; i++)                                     // Record the alignment infos and reference name for both reads
        {
            samOut[i]->hsp = (!i ? hitFirst->hit_hsps : hitSec->hit_hsps);
            samOut[i]->rname = safeStrdup((!i ? hitFirst->hit_def : hitSec->hit_def));
        }
        hitSec->hit_hsps = hitSec->hit_hsps->next;                  // Go to the second read's next HSP
        if (hitSec->hit_hsps != NULL)
            return hitRecord(hitFirst, hitSec, itSam, rVar);        // Record the other HitSec HSPs if there are any
        else
        { // If all the hitSec HSPs have been recorded with the current hitFirst HSP, go to the next hitFirst HSP
            hitFirst->hit_hsps = hitFirst->hit_hsps->next;
            if (hitFirst->hit_hsps != NULL)
            {
                hitSec->hit_hsps = rVar->hspTmp;                    // Go back to the first HSP of the second read
                rVar->hspTmp = NULL;
                itSam->samHits[countHit-1]->countHSPsec = 0;
                return hitRecord(hitFirst, hitSec, itSam, rVar);    // Record the other HitFirst HSPs if there are any
            }
        }
    }

    // If there are more reference hits for the first read
    if (hitFirst->next != NULL)
    {
        itSam->countHit++;
        itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
        itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
        if (hitSec == NULL)             // Single end or the end of the hitSec list have been reached
        {
            if (rVar->hitTmp == NULL)       // The second read is unmapped
                return hitRecord(hitFirst->next, NULL, itSam, rVar);
            else                            // The second read is mapped somewhere. hitTmp is used to go back to the first reference hit of the second read
                return hitRecord(hitFirst->next, rVar->hitTmp, itSam, rVar);
        }
        else                            // Paired end
            return hitRecord(hitFirst->next, hitSec->next, itSam, rVar);
    }

    // If the second read is mapped
    if (rVar->hitTmp != NULL)
    {
        if (!rVar->end)
        {
            hitSec = rVar->hitTmp;      // Go back to the first reference on which the second read is mapped
            rVar->end = 1;              // Mark that all the hits for the first read have been recorded
        }

        // If the second read has more than one reference hit
        if (hitSec->next != NULL)
        {
            // Go through all the previous hit records and stop if there is a reference name that matches the next reference
            for (rVar->tmpHitNb = 0; rVar->tmpHitNb < countHit; rVar->tmpHitNb++)
                if (!strcmp(hitSec->next->hit_def, itSam->samHits[rVar->tmpHitNb]->rsSam[0]->samOut[1]->rname)) break;

            // Occurs only if no match have been found
            if (rVar->tmpHitNb == countHit)
            {
                // Create a new Hit record
                itSam->countHit++;
                itSam->samHits = (SamHitPtr*) safeRealloc(itSam->samHits, itSam->countHit * sizeof(SamHitPtr)); // Append the Hit list
                itSam->samHits[itSam->countHit-1] = (SamHitPtr) safeCalloc(1, sizeof(SamHit));
            }
            return hitRecord(hitFirst, hitSec->next, itSam, rVar);
        }
    }
    return itSam;   // Return the whole superstructure when everything have been recorded about the read(s)
}


/************************************************************************************/
/*  Link the read's infos to its Blast results and create a SAM file                */
/************************************************************************************/
static IterationSamPtr iterationRecord(xmlTextReaderPtr reader, gzFile fp, gzFile fp2, AppParamPtr app)
{
    IterationSamPtr itSam = NULL;
    IterationPtr itFirst = NULL;
    IterationPtr itSec = NULL;
    RVarPtr rVar = (RVarPtr) safeCalloc(1, sizeof(RVar));   // Variables for hitRecord()

    // Paired end
    if (fp2 != NULL || app->inter)
    {
        itFirst = parseIteration(reader);       // Get the Blast results of the first in pair
        itSec = parseIteration(reader);         // Get the Blast results of the second in pair
        rVar->hitTmp = itSec->iteration_hits;   // Pointer to the first reference hit of the second in pair

        rVar->reads[0] = shortReadNext(fp);     // Get the first in pair infos from the first fastQ
        if (app->inter)                         // Get the second in pair infos either from the same fastQ if interleaved or from the second fastQ
            rVar->reads[1] = shortReadNext(fp);
        else
            rVar->reads[1] = shortReadNext(fp2);

        itSam = hitRecord(itFirst->iteration_hits, itSec->iteration_hits, itSam, rVar); // Record the results of the pair of reads together
        recordAnalysis(itSam, app);                                                     // Analyse the records
        printSam(itSam, app);                                                           // Print the two reads in the alignment section of the SAM file

        deallocIteration(itFirst);
        deallocIteration(itSec);
        shortReadFree(rVar->reads[0]);
        shortReadFree(rVar->reads[1]);
    }

    // Single end
    else
    {
        itFirst = parseIteration(reader);                               // Get the Blast results
        rVar->reads[0] = shortReadNext(fp);                             // Get the read's infos from the fastQ
        itSam = hitRecord(itFirst->iteration_hits, NULL, itSam, rVar);  // Record the results
        recordAnalysis(itSam, app);                                     // Analyse the records
        printSam(itSam, app);                                           // Print the read in the alignment section of the SAM file
        deallocIteration(itFirst);                                      // Dealloc the record structure
        shortReadFree(rVar->reads[0]);                                  // Dealloc the structure containing the read's infos
    }

    free(rVar);
    return itSam;
}


/************************************************************************************/
/*  Deallocation of the record structures                                           */
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
/*  Read group ID function                                                          */
/************************************************************************************/
// Extract the ID of the read group. Useful for Sam metadata
static char* readGroupID(char* readGroup)
{
    char* rgID = NULL;
    char* str = strstr(readGroup, "\tID:");
    str += 4;   // Set the pointer to the first char after the colon
    size_t i, str_length = strlen(str);
    int c = 0;

    // Count the number of char in the ID
    for (i = 0; i <= str_length && c != '\t'; i++)
        c = str[i];

    rgID = safeMalloc(i);
    memcpy(rgID, str, i-1);
    rgID[i-1] = '\0';
    return rgID;
}


/************************************************************************************/
/*  Main function of BlastSam.c                                                     */
/************************************************************************************/
int blastToSam(AppParamPtr app)
{
    xmlTextReaderPtr reader;
    gzFile fp = NULL, fp2 = NULL;
    BlastOutputPtr blastOP = NULL;
    IterationSamPtr itSam = NULL;

    reader = safeXmlNewTextReaderFilename(app->blastOut);       // Initialize the XML parsing

    safeXmlTextReaderRead(reader);                              // Get to the first line of the output

    if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
        ERROR("The document is not a Blast output\n", 1);

    blastOP = parseBlastOutput(reader);                         // Parse the blast output until it reaches the first iteration of results

    if (samHead(app) == 1)                                      // Print the Sam header
        ERROR("Error while printing the Sam header\n", 1);

    fp = safeGzOpen(app->fastq1, "r");                          // Open the first fastQ

    if (app->fastq2 != NULL)
        fp2 = safeGzOpen(app->fastq2, "r");                     // Open the second fastQ if there is one

    if (app->readGroup != NULL)
        app->readGroupID = readGroupID(app->readGroup);         // Extract the ID from the read group

    // Go through all the reads (iterations) of the XML output
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
