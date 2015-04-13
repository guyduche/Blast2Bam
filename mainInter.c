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

int main(int argc, char **argv)
{
	xmlTextReaderPtr reader;
	BlastOutputPtr blastOP = NULL;
	int evt = 1;
	int *cigar = NULL;
	int sizeCStr = 0;
	int i = 0;
	Iteration *itFor = NULL;
	Iteration *itRev = NULL;
	Hsp *hspFor = NULL;
	Hsp *hspRev = NULL;
	gzFile fp;
	kseq_t *seqFor;
	kseq_t *seqRev;

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

	seqFor = kseq_init(fp);
	seqRev = kseq_init(fp);

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{	
		itFor = parseIteration(reader);
		kseq_read(seqFor);
		itRev = parseIteration(reader);
		kseq_read(seqRev);

		if (!strcmp(seqFor->name.s, itFor->iteration_query_def) && !strcmp(seqRev->name.s, itRev->iteration_query_def))
		{
			if (itFor->iteration_hits->hit_hsps != NULL && itRev->iteration_hits->hit_hsps != NULL)
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
							printRead(itFor, hspFor, hspRev->hsp_hit_from, seqFor, sizeCStr);
							free(cigar);

							cigar = cigarStrBuilding(hspRev, itRev->iteration_query_len, &sizeCStr);
							printRead(itRev, hspRev, hspFor->hsp_hit_from, seqRev, sizeCStr);
							free(cigar);
						}
						hspRev = hspRev->next;
					}
					hspFor = hspFor->next;
				}

			else if ( // to be continued
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


void printRead(Iteration *it, Hsp *hsp, int posMate, kseq_t *seq, int sizeCStr)
{
	fprintf(stdout, "%s\tflag\t%s\t%d\t%d\t", it->iteration_query_def, it->iteration_hits->hit_def, hsp->hsp_hit_from, hsp->hsp_score);
	for (i = 0; i < sizeCStr; i++)
	{
		if (i % 2 == 0)
			fprintf(stdout, "%d", cigar[i]);
		else
			fprintf(stdout, "%c", cigar[i]);
	}	
	fprintf(stdout, "\t=\t%d\t%d\t%s\t%s\n", posMate, (hsp->hsp_hit_to)-(hsp->hsp_hit_from), seq->seq.s, seq->qual.s);
}









