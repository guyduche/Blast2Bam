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
#include "utils.h"
#include "blastSam.h"

/************************************************************************************/
/*  SAM header                                                                      */
/************************************************************************************/
// Extract the reference name
static char* refName(FILE* reader)
{
    int c = 0;
    size_t sizeStr = 0;
    char* name = NULL;

    c = fgetc(reader);                                          // Skip the ">"
    if (c != '>') return NULL;

    c = fgetc(reader);

    while (c != ' ' && c != '\t' && c != '\n' && c != EOF)      // Keep only the reference name
    {
        sizeStr++;
        name = (char*) safeRealloc(name, sizeStr+1);
        name[sizeStr-1] = (char) c;
        c = fgetc(reader);
    }

    name[sizeStr] = '\0';
    do c = fgetc(reader); while (c != '\n' && c != EOF);        // Skip the rest of the line

    return name;
}

// Extract the reference length
static int refLen(FILE* reader)
{
    int c = 0, len = 0, end = 0;

    while (!end && c != EOF)                                    // Until it is the end of the sequence and/or of the file
    {
        c = fgetc(reader);
        if (c == '\n')
        {
            c = fgetc(reader);
            if (c == '>')                                       // If there is another reference in the file
            {
                end = 1;
                fseek(reader, -1, SEEK_CUR);
            }
            else if (c == EOF || c == '\n');                    // If this is the end of the file or if there is another '\n'
            else len++;
        }
        else if (c == EOF);
        else len++;
    }
    return len;
}

// Print the SQ lines
static int sq_line(AppParamPtr app)
{
    FILE* reader;
    char* name = NULL;
    int len = 0;

    reader = safeFOpen(app->db, "r");                           // Open the reference file

    while (!feof(reader))
    {
        if ((name = refName(reader)) == NULL) return 1;         // Get the reference name
        len = refLen(reader);                                   // Get the length of the reference sequence

        fprintf(app->out, "@SQ\tSN:%s\tLN:%d\n", name, len);    // Print the SQ line in the SAM file
        free(name);
    }

    fclose(reader);
    return 0;
}

// Interpret the read group line: replace \t with '\t'
static char* interpretRG(AppParamPtr app)
{
    char *p, *q;
    char* str = safeStrdup(app->readGroup);

    for (p = q = str; *p != '\0'; p++, q++)
    {
        if (*p == '\\')
        {
            p++;
            if (*p == 't') *q = '\t';
        }
        else *q = *p;
    }
    *q = '\0';
    return str;
}

// Print the read group line
static int rg_line(AppParamPtr app)
{
    if (strstr(app->readGroup, "@RG") != app->readGroup) return 1;  // Verify that the read group begins with @RG
    if (strstr(app->readGroup, "\tID:") == NULL) return 1;          // Verify that the read group has an ID
    fprintf(app->out, "%s\n", app->readGroup);                      // Print the read group in the SAM header
    return 0;
}

// Print SAM header section
int samHead(AppParamPtr app)
{
    if (sq_line(app) == 1) return 1;                                // Print the SQ Line
    if (app->readGroup != NULL)                                     // If the read group is available and is correctly formatted, prints it
    {
        app->readGroup = interpretRG(app);                          // Interpret the read group
        if (rg_line(app) == 1) return 1;                            // If the read group is incorrectly formatted, returns an error
    }
    fprintf(app->out, "%s\n", app->pg_line);                        // Print the PG line
    return 0;
}


/************************************************************************************/
/*  Record analysis                                                                 */
/************************************************************************************/
// Function to filter the results by the alignment size
static int allowedToPrint(SamOutputPtr* samOut, int minLen, int countRec, int countRecGlobal, int countHSPsec, int countUnprintGlobal, int countUnprint, int countUnprintSec)
{
    int minLenAuto[2] = {0, 0};

    if (minLen != -1)                                                                                    // Option W is set
        minLenAuto[0] = minLenAuto[1] = minLen;
    else                                                                                                 // If not given, the minimum alignment length is calculated based on the read length
    {
        minLenAuto[0] = (samOut[0]->query->read_len / 5 > 20 ? samOut[0]->query->read_len / 5 : 20);        // Default filter: 20% of the read length must be aligned, with a minimum of 20
        if (samOut[1] != NULL)
            minLenAuto[1] = (samOut[1]->query->read_len / 5 > 20 ? samOut[1]->query->read_len / 5 : 20);
    }

    if (samOut[0]->hsp != NULL && samOut[0]->hsp->hsp_align_len > minLenAuto[0])                         // First read is mapped and alignment length is greater than the minimum length
    {
        if (samOut[1] != NULL && samOut[1]->hsp != NULL && samOut[1]->hsp->hsp_align_len < minLenAuto[1])   // Second read is mapped and alignment length is less than the minimum length
        {
            if (countUnprintSec != countHSPsec - 1)                                                             // Unless all the second read HSPs corresponding to the current first read HSP have been filtered, the record won't be printed
                return 0;
            else                                                                                                // If all the others have been filtered, print the record with second read considered unmapped
                samOut[1]->hsp = NULL;
        }
    }

    else                                                                                                 // First read is unmapped or alignment length less than the minimum length
    {
        if (samOut[1] != NULL && samOut[1]->hsp != NULL && samOut[1]->hsp->hsp_align_len > minLenAuto[1])   // Second read is mapped and alignment length greater than the minimum length
        {
            if (countUnprint < (countRec - countHSPsec))                                                        // If not last first read HSP, the record is not printed
                return 0;
            else                                                                                                // If last first read HSP or if first read unmapped, print with first read considered unmapped
                samOut[0]->hsp = NULL;
        }
        else                                                                                                // Single end or second read unmapped or alignment length less than the minimum length
        {
            if (countUnprintGlobal != countRecGlobal -1)                                                        // If not last record, it is not printed
                return 0;
            else                                                                                                // If last record, print with both reads considered unmapped
            {
                samOut[0]->hsp = NULL;
                if (samOut[1] != NULL)
                    samOut[1]->hsp = NULL;
            }
        }
    }
    return 1;
}

// Analyse the records to flag the best and flag the ones unfit to be printed in the SAM file
void recordAnalysis(IterationSamPtr itSam, AppParamPtr app)
{
    int i = 0, j = 0, len0 = 0, len1 = 0, countUnprint = 0, countUnprintGlobal = 0, countUnprintSec = 0;
    double score = 0.0, bestScore = 10.0;
    RecordSamPtr bestRecord = NULL;
    RecordSamPtr curRecord = NULL;

    for (i = 0; i < itSam->countHit; i++)               // Go through all the reference hits
    {
        for (j = 0; j < itSam->samHits[i]->countRec; j++)   // Go through all the records
        {
            curRecord = itSam->samHits[i]->rsSam[j];
            if (!i && !j) bestRecord = curRecord;
            score = 0.0;

            // CountUnprintSec is reset for each new HSP first
            if (itSam->samHits[i]->countHSPsec != 0 && j % itSam->samHits[i]->countHSPsec == 0)
                countUnprintSec = 0;

            // Filter the records
            if (!allowedToPrint(curRecord->samOut, app->minLen, itSam->samHits[i]->countRec, itSam->countRecGlobal, itSam->samHits[i]->countHSPsec, countUnprintGlobal, countUnprint, countUnprintSec))
            {
                curRecord->doNotPrint = 1;
                countUnprint++;
                countUnprintSec++;
                countUnprintGlobal++;
                continue;
            }

            // The first in pair is mapped
            if (curRecord->samOut[0]->hsp != NULL)
            {
                score += curRecord->samOut[0]->hsp->hsp_evalue;
                // The second in pair is mapped
                if (curRecord->samOut[1] != NULL && curRecord->samOut[1]->hsp != NULL)
                {
                    score += curRecord->samOut[1]->hsp->hsp_evalue;
                    curRecord->tlen = abs(curRecord->samOut[0]->hsp->hsp_hit_from - curRecord->samOut[1]->hsp->hsp_hit_from) + 1;       // Tlen
                    len0 = abs(curRecord->samOut[0]->hsp->hsp_hit_to - curRecord->samOut[0]->hsp->hsp_hit_from) + 1;                    // Length of the first in pair alignment on the reference
                    len1 = abs(curRecord->samOut[1]->hsp->hsp_hit_to - curRecord->samOut[1]->hsp->hsp_hit_from) + 1;                    // Length of the second in pair alignment on the reference
                    if (curRecord->tlen > 3 * (len0 >= len1 ? len0 : len1))                                                             // 3 times the alignment length of the longest read
                    {                                                                                                                   // Pair of reads not properly aligned
                        if (countUnprintGlobal != (itSam->countRecGlobal -1))                                                           // If not the last record, go to the next record
                        {
                            curRecord->doNotPrint = 1;
                            countUnprint++;
                            countUnprintSec++;
                            countUnprintGlobal++;
                            continue;
                        }
                        else                                                                                                            // If last record, both reads are considered unmapped
                            curRecord->samOut[0]->hsp = curRecord->samOut[1]->hsp = NULL;
                    }
                }
                // The second in pair is unmapped
                else
                {
                    score += 0.5;
                    curRecord->tlen = 0;
                }
            }
            // The first in pair is unmapped
            else
            {
                score += 0.5;
                curRecord->tlen = 0;
                // The second in pair is mapped
                if (curRecord->samOut[1] != NULL && curRecord->samOut[1]->hsp != NULL)
                    score += curRecord->samOut[1]->hsp->hsp_evalue;
                else
                    score += 0.5;
            }

            if (score < bestScore)
            {
                bestScore = score;
                bestRecord = curRecord;
            }
        }
        countUnprint = 0;
    }
    bestRecord->best = 1;
}


/************************************************************************************/
/*  Build the CIGAR string corresponding to the read alignment                      */
/************************************************************************************/
// Put the count and the symbol in the CIGAR structure
#define CSTRMACRO(TAG, FUNCT) { \
        COUNTCHARMACRO(FUNCT); \
        cigar->elements[cigar->size-1]->count = count; \
        cigar->elements[cigar->size-1]->symbol = TAG;}

// Count the number of alignment event
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

    if (hsp->hsp_query_from > 1)                // 5' Soft clipping
    {
        cigar->size++;
        cigar->elements = (CigarElementPtr*) safeCalloc(cigar->size, sizeof(CigarElementPtr));  // Create an array of elements for the CIGAR string
        cigar->elements[0] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));             // Create the first element of the array
        cigar->elements[0]->count = hsp->hsp_query_from - 1;                                    // Count = first alignment position in the query sequence - 1
        cigar->elements[0]->symbol = 'S';
    }

    for (pos = 0; pos < (hsp->hsp_align_len); pos++)
    {
        cigar->size++;
        cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr));   // Create or append the array of elements of the CIGAR string
        cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));                     // Create a new element

        if (hsp->hsp_hseq[pos] == '-')
            CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-'))                             // Count the number of insertions

        else if (hsp->hsp_qseq[pos] == '-')
            CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-'))                             // Count the number of deletions

        else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
            CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos]))              // Count the number of matches

        else
            CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos]))    // Count the number of mismatches

        if (cigar->elements[cigar->size-1]->count >= 100 && cigar->elements[cigar->size-1]->symbol == 'D')
            cigar->elements[cigar->size-1]->symbol = 'N';                           // If there is more than a hundred deletion at a time, it is considered a skipped region

        if (cigar->elements[cigar->size-1]->symbol != '=')
            cigar->nbDiff += cigar->elements[cigar->size-1]->count;                 // Count the number of I, D/N and X for the NM tag
    }

    if ((queryLength - hsp->hsp_query_to) > 0)  // 3' Soft clipping
    {
        cigar->size++;
        cigar->elements = (CigarElementPtr*) safeRealloc(cigar->elements, cigar->size * sizeof(CigarElementPtr));
        cigar->elements[cigar->size-1] = (CigarElementPtr) safeCalloc(1, sizeof(CigarElement));
        cigar->elements[cigar->size-1]->count = queryLength - hsp->hsp_query_to;    // Count = length of the read - position of the last aligned base on the query sequence
        cigar->elements[cigar->size-1]->symbol = 'S';
    }

    return cigar;
}


/************************************************************************************/
/*  Get the reverse complement of a sequence                                        */
/************************************************************************************/
static char* revStr(char* oldStr)
{
    size_t i = 0;
    int l = strlen(oldStr);
    char* newStr = safeCalloc(l+1, sizeof(char));

    for (--l; l >= 0; l--, i++)                     // Read oldStr in reverse and newStr straight
    {
        switch (oldStr[l])                          // For each char in oldStr, put the complement in the newStr
        {
            case 'A': case 'a' : newStr[i] = 'T'; break;
            case 'T': case 't' : newStr[i] = 'A'; break;
            case 'C': case 'c' : newStr[i] = 'G'; break;
            case 'G': case 'g' : newStr[i] = 'C'; break;
            default: newStr[i] = oldStr[l]; break;  // If not in {ATCG}, keep the same char
        }
    }
    return newStr;
}


/************************************************************************************/
/*  Extract the beginning position in the reference name, if there is one           */
/************************************************************************************/
// Search for the first colon in the reference name and extract the number directly after
static int firstPosRef(const char* rname)
{
    char* colon = strchr(rname, ':');
    if (colon == NULL)
        return 0;
    return strtol(colon + 1, NULL, 10);
}


/************************************************************************************/
/*  Print SAM alignment section                                                     */
/************************************************************************************/
// SAM flags
#define SAM_PAIRED 0x1              // Paired end
#define SAM_PROPER_PAIR 0x2         // Read mapped in a proper pair
#define SAM_UNMAP 0x4               // Read unmapped
#define SAM_MUNMAP 0x8              // Mate unmapped
#define SAM_REVERSE 0x10            // Read mapped on the reverse strand
#define SAM_MREVERSE 0x20           // Mate mapped on the reverse strand
#define SAM_READF 0x40              // Read is first in pair
#define SAM_READL 0x80              // Read is last in pair
#define SAM_SECONDARY 0x100         // Not primary alignment
#define SAM_QCFAIL 0x200            // Failed control quality
#define SAM_DUP 0x400               // Optical or PCR duplicate
#define SAM_SUPPLEMENTARY 0x800     // Supplementary alignment

// Print the CIGAR str straight or reverse depending on the flag
static void printCigarStr (AppParamPtr app, CigarElementPtr* cigElements, size_t size, int flag)
{
    int i = 0;
    if (flag & SAM_REVERSE)                 // If the sequence is mapped on the reverse strand
        for (i = size - 1; i >= 0; i--)     // Go through the CIGAR elements in reverse order
        {
            fprintf(app->out, "%d%c", cigElements[i]->count, cigElements[i]->symbol);
            free(cigElements[i]);
        }

    else                                    // If the sequence is mapped on the forward strand
        for (i = 0; i < size; i++)          // Go through the CIGAR element in straight order
        {
            fprintf(app->out, "%d%c", cigElements[i]->count, cigElements[i]->symbol);
            free(cigElements[i]);
        }

    free(cigElements);                      // The CIGAR string is freed just after being printed
}

// Structure that contains the infos of one line of SAM alignment section
typedef struct SamLine
{
    char* readName;         // Col 1    QNAME
    char* refName;          // Col 3    RNAME
    char* rnext;            // Col 7    RNEXT
    char* seq;              // Col 10   SEQ
    char* qual;             // Col 11   QUAL
    int flag;               // Col 2    FLAG
    int pos;                // Col 4    POS
    int mapq;               // Col 5    MAPQ
    int pnext;              // Col 8    PNEXT
    int tlen;               // Col 9    TLEN
    int blastScore;         // Col 12   Metadata AS
    double blastBitScore;   // Col 12   Metadata XB
    double blastEvalue;     // Col 12   Metadata XE
    CigarPtr cigarStr;      // Col 6    CIGAR
} SamLine, *SamLinePtr;

// Print one line of the SAM alignment section
static void printSamLine (AppParamPtr app, SamLinePtr samLine)
{
    int i = 0;
    size_t qualLen = 0;
    char* seq = NULL;

    fputs(samLine->readName, app->out);                                                             // Print QNAME
    fprintf(app->out, "\t%d\t", samLine->flag);                                                     // Print FLAG
    fputs(samLine->refName, app->out);                                                              // Print RNAME
    fprintf(app->out, "\t%d\t%d\t", samLine->pos, samLine->mapq);                                   // Print POS and MAPQ

    if (samLine->cigarStr != NULL)
    {
        printCigarStr (app, samLine->cigarStr->elements, samLine->cigarStr->size, samLine->flag);   // Print CIGAR
        free(samLine->cigarStr);
    }
    else
        fputc('*', app->out);

    fputc('\t', app->out);
    fputs(samLine->rnext, app->out);                                                                // Print RNEXT
    fprintf(app->out, "\t%d\t%d\t", samLine->pnext, samLine->tlen);                                 // Print PNEXT and TLEN

    if (samLine->flag & SAM_REVERSE)                                                                // If mapped on the reverse strand, print the reverse complement SEQ and the reverse of QUAL
    {
        seq = revStr(samLine->seq);
        fputs(seq, app->out);
        fputc('\t', app->out);

        qualLen = strlen(samLine->qual);
        for (i = qualLen - 1; i >= 0; i--)
            fputc(samLine->qual[i], app->out);
    }
    else                                                                                            // If mapped on the forward strand, print SEQ and QUAL as they came out of the sequencer
    {
        fputs(samLine->seq, app->out);
        fputc('\t', app->out);
        fputs(samLine->qual, app->out);
    }
                                                                                                    // Metadata
    if (samLine->cigarStr != NULL)
        fprintf(app->out, "\tNM:i:%d", samLine->cigarStr->nbDiff);                                  // Print NM tag

    if (app->readGroupID != NULL)
        fprintf(app->out, "\tRG:Z:%s", app->readGroupID);                                           // Print RG tag

    fprintf(app->out, "\tAS:i:%d", samLine->blastScore);                                            // Print AS tag
    fprintf(app->out, "\tXB:f:%g", samLine->blastBitScore);                                         // Print XB tag
    fprintf(app->out, "\tXE:Z:%.3g", samLine->blastEvalue);                                         // Print XE tag. Z and not f because BAM files can't handle those numbers
    fputc('\n', app->out);
}

// Print the alignment section
void printSam(IterationSamPtr itSam, AppParamPtr app)
{
    int i = 0, j = 0, k = 0, invk = 0, posRef = 0;
    SamOutputPtr samOut[2] = {NULL, NULL};
    SamLinePtr samLine = NULL;

    for (i = 0; i < itSam->countHit; i++)               // Go through all the reference hits
    {
        for (j = 0; j < itSam->samHits[i]->countRec; j++)   // Go through all the records
        {
            if (itSam->samHits[i]->rsSam[j]->doNotPrint) continue;

            for (k = 0; k < 2; k++)                             // Print the first in pair and then the second
            {
                invk = (k ? 0 : 1);
                samOut[k] = itSam->samHits[i]->rsSam[j]->samOut[k];
                samOut[invk] = itSam->samHits[i]->rsSam[j]->samOut[invk];

                if (samOut[k] == NULL) continue;                        // When on the second in pair part of the record, if single end, go to the next record

                if (!k)
                {
                    if (app->posOnChr && samOut[k]->rname != NULL) posRef = firstPosRef(samOut[k]->rname);  // Extract the position of the reference from its name (-z)
                    else posRef = 0;
                }

                samLine = (SamLinePtr) safeCalloc(1, sizeof(SamLine));  // Create a new line for the alignment section of the SAM file

                // Paired end
                if (samOut[invk] != NULL)
                {
                    samLine->flag |= SAM_PAIRED | (!k ? SAM_READF : SAM_READL);

                    // The mate is mapped
                    if (samOut[invk]->hsp != NULL)
                    {
                        samLine->flag |= (samOut[invk]->hsp->hsp_hit_to < samOut[invk]->hsp->hsp_hit_from ? SAM_MREVERSE : 0);              // Mate mapped on the reverse strand ?
                        samLine->pnext = (samLine->flag & SAM_MREVERSE ? samOut[invk]->hsp->hsp_hit_to : samOut[invk]->hsp->hsp_hit_from);  // PNEXT is the leftmost position of the mate alignment on the reference
                        samLine->pnext += posRef;                                                                                           // Adjust the position to the first position of the reference (-z)

                        // The read is mapped
                        if (samOut[k]->hsp != NULL)
                        {
                            samLine->rnext = "=";
                            samLine->flag |= SAM_PROPER_PAIR;
                        }

                        // The read is unmapped
                        else
                        {
                            samLine->refName = shortName(samOut[invk]->rname);                                                              // According to SAM specs, unmapped reads should have
                            samLine->pos = samLine->pnext;                                                                                  // the RNAME and POS of their mate
                            samLine->rnext = "=";
                        }
                    }

                    // The mate is unmapped
                    else
                    {
                        samLine->flag |= SAM_MUNMAP;
                        samLine->rnext = "*";
                    }
                }

                // Single end
                else
                    samLine->rnext = "*";

                // Put the read infos in samLine
                samLine->readName = samOut[k]->query->name;
                samLine->seq = samOut[k]->query->seq;
                samLine->qual = samOut[k]->query->qual;
                if (samLine->qual == NULL)
                    samLine->qual = "*";

                // The read is mapped
                if (samOut[k]->hsp != NULL)
                {
                    if (!itSam->samHits[i]->rsSam[j]->best)                                                                     // This is not the best record of this read
                    {
                        samLine->flag |= SAM_SECONDARY;                                                                         // Secondary alignment
                        samLine->seq = "*";
                        samLine->qual = "*";                                                                                    // Reduce file size
                    }
                    samLine->flag |= (samOut[k]->hsp->hsp_hit_to < samOut[k]->hsp->hsp_hit_from ? SAM_REVERSE : 0);             // The read is mapped on the reverse strand ?
                    samLine->refName = shortName(samOut[k]->rname);                                                             // Put a short version of the reference name in samLine
                    samLine->pos = (samLine->flag & SAM_REVERSE ? samOut[k]->hsp->hsp_hit_to : samOut[k]->hsp->hsp_hit_from);   // POS is the leftmost position of the read alignment on the reference
                    samLine->pos += posRef;                                                                                     // Adjust the position to the first position of the reference (-z)
                    samLine->cigarStr = cigarStrBuilding(samOut[k]);                                                            // Build the CIGAR string
                    samLine->mapq = 60;                                                                                         // MAPQ
                    samLine->tlen = itSam->samHits[i]->rsSam[j]->tlen;                                                          // TLEN: distance between the first alignment position of both reads
                    if (samLine->tlen && samLine->flag & SAM_REVERSE)
                        samLine->tlen = -(samLine->tlen);
                    samLine->blastScore = samOut[k]->hsp->hsp_score;                                                            // AS metadata tag: Blast score
                    samLine->blastBitScore = samOut[k]->hsp->hsp_bit_score;                                                     // XB metadata tag: Blast bit score
                    samLine->blastEvalue = samOut[k]->hsp->hsp_evalue;                                                          // XE metadata tag: Blast E-Value

                    // The mate is unmapped
                    if (samOut[invk] != NULL && samOut[invk]->hsp == NULL)
                    {
                        samLine->rnext = "=";
                        samLine->pnext = samLine->pos;
                    }
                }

                // The read is unmapped
                else
                {
                    samLine->flag |= SAM_UNMAP;
                    if (samOut[invk] != NULL && samOut[invk]->hsp != NULL && !itSam->samHits[i]->rsSam[j]->best)                // Read unmapped, mate mapped but secondary: read flagged secondary
                    {
                        samLine->flag |= SAM_SECONDARY;
                        samLine->seq = "*";
                        samLine->qual = "*";                                                                                    // Reduce file size
                    }
                    if (samLine->refName == NULL)
                        samLine->refName = "*";
                }

                // Print a new line in SAM alignment section
                printSamLine (app, samLine);
                free(samLine);
            }
        }
    }
}
