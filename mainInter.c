#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "cigar.h"
#include "macro.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct Sequences_t
{
	kseq_t *seqFor;
	kseq_t *seqRev;
} Sequences;

int main(int argc, char **argv)
{
	xmlTextReaderPtr reader;
	BlastOutputPtr blastOP = NULL;
	int evt = 1;
	int *cigar = NULL;
	int sizeCStr = 0;
	int flag = 0;
	Iteration *itFor = NULL;
	Iteration *itRev = NULL;
	Hsp *hspFor = NULL;
	Hsp *hspRev = NULL;
	Hsp *hspCur = NULL;
	gzFile fp;
	kseq_t *seq;
	Sequences *sequences = NULL;

	if (argc < 3)
		ERROR("Wrong number of arguments\n", NULL, EXIT_FAILURE)
	
	reader = xmlNewTextReaderFilename(argv[1]);
	
	if (reader == NULL)
		ERROR("Unable to open the XML file\n", NULL, EXIT_FAILURE)
	
	evt = xmlTextReaderRead(reader);	

	if (evt == -1)
		ERROR("Error while reading the first node\n", xmlFreeTextReader(reader), EXIT_FAILURE)

	if (xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "BlastOutput"))
		ERROR("The document is not a Blast output\n", xmlFreeTextReader(reader), EXIT_FAILURE)

	evt = xmlTextReaderRead(reader);

	blastOP = parseBlastOutput(reader);

	fp = gzopen(argv[2], "r");

	if (fp == NULL)
		ERROR("Unable to open the FastQ\n", xmlFreeTextReader(reader), EXIT_FAILURE)

	seq = kseq_init(fp);
	TABBUILDINGMACRO(sequences, 1, Sequences, EXIT_FAILURE)
	TABBUILDINGMACRO(sequences->seqFor, 1, kseq_t, EXIT_FAILURE)
	TABBUILDINGMACRO(sequences->seqRev, 1, kseq_t, EXIT_FAILURE)

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{	
		itFor = parseIteration(reader);
		kseq_read(seq);
		// problème de l'interleaved à résoudre !!!!!
		itRev = parseIteration(reader);
		kseq_read(seqRev);

		if (!strcmp(seqFor->name.s, itFor->iteration_query_def) && !strcmp(seqRev->name.s, itRev->iteration_query_def))
		{
			if (itFor->iteration_hits->hit_hsps == NULL && itRev->iteration_hits->hit_hsps == NULL)
			{
				printRead(itFor->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqFor->seq.s, seqFor->qual.s);
				printRead(itRev->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seqRev->seq.s, seqRev->qual.s);
			}

			else if (itFor->iteration_hits->hit_hsps == NULL)
			{
				hspFor = itFor->iteration_hits->hit_hsps;
				hspCur = hspFor;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspFor->hsp_score)
						hspFor = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(hspFor, itFor->iteration_query_len, &sizeCStr);
				printRead(itFor->iteration_query_def, flag, itFor->iteration_hits->hit_def, hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", 0, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq.s, seqFor->qual.s);
				free(cigar);

				printRead(itRev->iteration_query_def, flag, itFor->iteration_hits->hit_def, 0, 0, NULL, 0, "=", hspFor->hsp_hit_from, 0, seqRev->seq.s, seqRev->qual.s);
			}

			else if (itRev->iteration_hits->hit_hsps == NULL)
			{
				hspRev = itRev->iteration_hits->hit_hsps;
				hspCur = hspRev;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hspRev->hsp_score)
						hspRev = hspCur;
					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(hspRev, itRev->iteration_query_len, &sizeCStr);
				printRead(itRev->iteration_query_def, flag, itRev->iteration_hits->hit_def, hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", 0, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq.s, seqRev->qual.s);
				free(cigar);

				printRead(itFor->iteration_query_def, flag, itRev->iteration_hits->hit_def, 0, 0, NULL, 0, "=", hspRev->hsp_hit_from, 0, seqFor->seq.s, seqFor->qual.s);
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
						if ((300 < abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to)) && (abs(hspRev->hsp_hit_to - hspFor->hsp_hit_to) < 600))
						{
							cigar = cigarStrBuilding(hspFor, itFor->iteration_query_len, &sizeCStr);
							printRead(itFor->iteration_query_def, flag, itFor->iteration_hits->hit_def, hspFor->hsp_hit_from, hspFor->hsp_score, cigar, sizeCStr, "=", hspRev->hsp_hit_from, (hspFor->hsp_hit_to)-(hspFor->hsp_hit_from), seqFor->seq.s, seqFor->qual.s);
							free(cigar);

							cigar = cigarStrBuilding(hspFor, itFor->iteration_query_len, &sizeCStr);
							printRead(itRev->iteration_query_def, flag, itRev->iteration_hits->hit_def, hspRev->hsp_hit_from, hspRev->hsp_score, cigar, sizeCStr, "=", hspFor->hsp_hit_from, (hspRev->hsp_hit_to)-(hspRev->hsp_hit_from), seqRev->seq.s, seqRev->qual.s);
							free(cigar);
						}

						hspRev = hspRev->next;
					}
					hspFor = hspFor->next;
				}
			}

			else if (itFor->iteration_hits->hit_hsps == NULL && itRev->iteration_hits->hit_hsps == NULL)
			{
				printRead(
fprintf(stdout, "%s\tflag\t*\t0\t0\t*\trnext\tpnext\t0\t%s\t%s\n", it->iteration_query_def, seq->seq.s, seq->qual.s);
			}
		
			else
				fprintf(stdout, "%s\tflag\t*\t0\t0\t*\trnext\tpnext\t0\t%s\t%s\n", it->iteration_query_def, seq->seq.s, seq->qual.s);
		}
		deallocIteration(it);
	}
	kseq_destroy(seq);
	gzclose(fp);
	deallocBlastOutput(blastOP);
	
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	
	return EXIT_SUCCESS;
}


void printRead(char *qname, int flag, char *rname, int pos, int mapq, char *cigar, int sizeCStr, char *rnext, int pnext, int tlen, char *seq, char *qual)
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



