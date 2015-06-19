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
#ifndef BLASTSAM_H_INCLUDED
#define BLASTSAM_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include "shortRead.h"
#include "parseXML.h"

/************************************************************************************/
/*	Program parameters																*/
/************************************************************************************/
typedef struct AppParam
{
	int inter;					// Interleaved option (-p)
	int posOnChr;				// Position based on the reference name (-z)
	int minLen;					// Minimum alignment length accepted (-W)
	char* pg_line;				// @PG line of SAM header
	char* readGroup;			// Read group (-R)
	char* readGroupID;			// Read group ID
	char* blastOut;				// Blast XML output
	char* db;					// Reference file
	char* fastq1;				// First fastQ
	char* fastq2;				// Second fastQ
	FILE* out;					// Output stream (-o): default is stdout
}AppParam, *AppParamPtr;


/************************************************************************************/
/*	CIGAR string structures															*/
/************************************************************************************/
// Sub-unit of the CIGAR string
typedef struct CigarElement
{
	int count;					// Number of times an alignment event is encountered before another one occurs
	char symbol;				// Comprised in {I, D, N, X, =}
} CigarElement,*CigarElementPtr;

// CIGAR string
typedef struct Cigar
{
	CigarElementPtr* elements;	// Array of CIGAR sub-units
	int nbDiff;					// Number of differences between the read and the reference, clipped bases excluded
	size_t size;				// Size of the CIGAR string
} Cigar, *CigarPtr;


/************************************************************************************/
/*	Record structures																*/
/************************************************************************************/
// Sub-unit of RecordSam : contains the infos concerning the read and its alignment
typedef struct SamOutput
{
	ShortReadPtr query;			// Read infos from the fastQ (name, seq and qual)
	HspPtr hsp;					// Alignment infos from Blast
	char* rname;				// Name of the reference on which the read is mapped
} SamOutput, *SamOutputPtr;

// Sub-unit of SamHit : contains a paired alignment record
typedef struct RecordSam
{
	SamOutputPtr samOut[2];		// Array containing the reads and their alignment infos. 0: first in pair; 1: second in pair
} RecordSam, *RecordSamPtr;

// Sub-unit of IterationSam : contains the alignments for a single reference
// Each possible alignment for the first in pair is recorded along each of the possible alignment for the second in pair
typedef struct SamHit
{
	RecordSamPtr* rsSam;		// Array containing the paired records
	size_t countRec;			// Number of paired records (= number of first in pair HSPs x number of second in pair HSPs)
	size_t countHSPsec;			// Number of HSP found for the second in pair (useful for option -W)
} SamHit, *SamHitPtr;

// Structure containing all the alignments found by Blast for a read (single end) or a pair of reads (paired end)
typedef struct IterationSam
{
	SamHitPtr* samHits;			// Array containing the hits found for each reference
	size_t countHit;			// Number of hits (number of references on which an alignment has been found for the read)
} IterationSam, *IterationSamPtr;



/************************************************************************************/
/*	Prototypes																		*/
/************************************************************************************/
int blastToSam(AppParamPtr app);
int samHead(AppParamPtr app);
void printSam(IterationSamPtr itSam, AppParamPtr app);


#endif // BLASTSAM_H_INCLUDED

