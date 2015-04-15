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
	Iteration *it = NULL;
	Hsp *hsp = NULL;
	Hsp *runningHsp = NULL;
	gzFile fp;
	kseq_t *seq;

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

	seq = kseq_init(fp);

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		it = parseIteration(reader);
		kseq_read(seq);
		if (!strcmp(seq->name.s, it->iteration_query_def))
		{
			if (it->iteration_hits->hit_hsps != NULL)
			{			
				hsp = it->iteration_hits->hit_hsps;
				runningHsp = it->iteration_hits->hit_hsps;

				while (runningHsp != NULL)
				{
					if (runningHsp->hsp_score > hsp->hsp_score)
						hsp = runningHsp;

					runningHsp = runningHsp->next;
				}

				cigar = cigarStrBuilding(cigar, hsp, it->iteration_query_len, &sizeCStr);
				fprintf(stdout, "%s\tflag\t%s\t%d\t%d\t", it->iteration_query_def, it->iteration_hits->hit_def, hsp->hsp_hit_from, hsp->hsp_score);
				for (i = 0; i < sizeCStr; i++)
				{
					if (i % 2 == 0)
						fprintf(stdout, "%d", cigar[i]);
					else
						fprintf(stdout, "%c", cigar[i]);
				}	
				fprintf(stdout, "\trnext\tpnext\t%d\t%s\t%s\n", (hsp->hsp_hit_to)-(hsp->hsp_hit_from), seq->seq.s, seq->qual.s);
			}
		
			else
				fprintf(stdout, "%s\tflag\t*\t0\t0\t*\trnext\tpnext\t0\t%s\t%s\n", it->iteration_query_def, seq->seq.s, seq->qual.s);
		}
		deallocIteration(it);
	}
	free(cigar);
	kseq_destroy(seq);
	gzclose(fp);
	deallocBlastOutput(blastOP);
	
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	
	return EXIT_SUCCESS;
}

