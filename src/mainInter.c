
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "cigar.h"
#include "macro.h"
#include "shortRead.h"


void printRead(char *qname, int flag, char *rname, int pos, int mapq, int *cigar, int sizeCStr, char *rnext, int pnext, int tlen, char *seq, char *qual)
{
	int i = 0;
	fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t", qname, flag, rname, pos, mapq);
	if (cigar != NULL)
		for (i = 0; i < sizeCStr; i++)
			fprintf(stdout,(i % 2 == 0 ? "%d":"%c"), cigar[i]);
	else
		fprintf(stdout, "*");
	fprintf(stdout, "\t%s\t%d\t%d\t%s\t%s\n", rnext, pnext, tlen, seq, qual);
}

int main(int argc, char **argv)
{
	xmlTextReaderPtr reader;
	BlastOutputPtr blastOP = NULL;
	int evt = 1;
	int *cigar = NULL;
	int sizeCStr = 0;
	int flag = 0;
	IterationPtr itFor = NULL;
	IterationPtr itRev = NULL;
	Hsp *hspFor = NULL;
	Hsp *hspRev = NULL;
	Hsp *hspCur = NULL;
	ShortReadPtr seqFor = NULL;
	ShortReadPtr seqRev = NULL;
	gzFile fp;
	int side = 0;
	ShortReadAndIteration shortReadsIter[2];

	if (argc < 3)
		ERROR("Wrong number of arguments\n", EXIT_FAILURE)
	
	reader = xmlNewTextReaderFilename(argv[1]);
	
	if (reader == NULL)
		ERROR("Unable to open the XML file\n", EXIT_FAILURE)
	
	evt = xmlTextReaderRead(reader);	

	if (evt == -1)
		ERROR("Error while reading the first node\n", EXIT_FAILURE)

	if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
		ERROR("The document is not a Blast output\n", EXIT_FAILURE)

	evt = xmlTextReaderRead(reader);

	blastOP = parseBlastOutput(reader);

	fp = gzopen(argv[2], "r");

	if (fp == NULL)
		ERROR("Unable to open the FastQ\n", EXIT_FAILURE)

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		for(side = 0; side < 2; ++side)
		{
			shortReadsIter[side].blast = parseIteration(reader);
			shortReadsIter[side].seq = ShortReadNext(fp);
		}

		itFor = shortReadsIter[0].blast;
		itRev = shortReadsIter[1].blast;
		seqFor = shortReadsIter[0].seq;
		seqRev = shortReadsIter[1].seq;

		if (!strcmp(seqFor->name, itFor->iteration_query_def) && !strcmp(seqRev->name, itRev->iteration_query_def))
		{
			if (itFor->iteration_hits->hit_hsps == NULL && itRev->iteration_hits->hit_hsps == NULL)
			{
				printRead(itFor->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqFor->seq, seqFor->qual);
				printRead(itRev->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqRev->seq, seqRev->qual);
			}

			else if (itRev->iteration_hits->hit_hsps == NULL)
			{
				hspFor = itFor->iteration_hits->hit_hsps;
				hspCur = hspFor;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspFor->hsp_score)
						hspFor = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hspFor, itFor->iteration_query_len, &sizeCStr);
				printRead(itFor->iteration_query_def, flag, itFor->iteration_hits->hit_def, hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", 0, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);

				printRead(itRev->iteration_query_def, flag, itFor->iteration_hits->hit_def, 0, 0, NULL, 0, "=", hspFor->hsp_hit_from, 0, seqRev->seq, seqRev->qual);
			}

			else if (itFor->iteration_hits->hit_hsps == NULL)
			{
				hspRev = itRev->iteration_hits->hit_hsps;
				hspCur = hspRev;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspRev->hsp_score)
						hspRev = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hspRev, itRev->iteration_query_len, &sizeCStr);
				printRead(itRev->iteration_query_def, flag, itRev->iteration_hits->hit_def, hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", 0, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq, seqRev->qual);

				printRead(itFor->iteration_query_def, flag, itRev->iteration_hits->hit_def, 0, 0, NULL, 0, "=", hspRev->hsp_hit_from, 0, seqFor->seq, seqFor->qual);
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
						if ((0 < abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to)) && (abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to) < 1000000))
						{
							cigar = cigarStrBuilding(cigar, hspFor, itFor->iteration_query_len, &sizeCStr);
							printRead(itFor->iteration_query_def, flag, itFor->iteration_hits->hit_def, hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", hspRev->hsp_hit_from, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq, seqFor->qual);

							cigar = cigarStrBuilding(cigar, hspRev, itRev->iteration_query_len, &sizeCStr);
							printRead(itRev->iteration_query_def, flag, itRev->iteration_hits->hit_def, hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", hspFor->hsp_hit_from, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq, seqRev->qual);
						}

						hspRev = hspRev->next;
					}
					hspFor = hspFor->next;
				}
			}
		}
		
		for(side = 0; side < 2; ++side)
		{
			deallocIteration(shortReadsIter[side].blast);
			ShortReadFree(shortReadsIter[side].seq);
		}
	}
	if (cigar != NULL)
		free(cigar);

	gzclose(fp);
	deallocBlastOutput(blastOP);
	
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	
	return EXIT_SUCCESS;
}





