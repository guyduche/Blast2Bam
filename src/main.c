
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
	IterationPtr it = NULL;
	Hsp *hsp = NULL;
	Hsp *hspCur = NULL;
	gzFile fp;
	ShortReadPtr seq = NULL;

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
		it = parseIteration(reader);
		seq = ShortReadNext(fp);
		if (!strcmp(seq->name, it->iteration_query_def))
		{
			if (it->iteration_hits->hit_hsps != NULL)
			{			
				hsp = it->iteration_hits->hit_hsps;
				hspCur = it->iteration_hits->hit_hsps;

				while (hspCur != NULL)
				{
					if (hspCur->hsp_score > hsp->hsp_score)
						hsp = hspCur;

					hspCur = hspCur->next;
				}

				cigar = cigarStrBuilding(cigar, hsp, it->iteration_query_len, &sizeCStr);
				printRead(it->iteration_query_def, flag, it->iteration_hits->hit_def, hsp->hsp_hit_from, hsp->hsp_score, cigar, sizeCStr, "*", 0, (hsp->hsp_hit_to)-(hsp->hsp_hit_from), seq->seq, seq->qual);
			}
		
			else
				printRead(it->iteration_query_def, 0, "*", 0, 0, NULL, 0, "*", 0, 0, seq->seq, seq->qual);
		}
		deallocIteration(it);
	}
	if (cigar != NULL)
		free(cigar);
	ShortReadFree(seq);
	gzclose(fp);
	deallocBlastOutput(blastOP);
	
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	
	return EXIT_SUCCESS;
}

