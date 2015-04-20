
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlreader.h>
#include <zlib.h>
#include "parseXML.h"
#include "blastSam.h"
#include "utils.h"
#include "shortRead.h"

//ALLER VOIR CODE JENNIFER get_longopt / getopt

int main(int argc, char** argv)
{
	xmlTextReaderPtr reader;
	BlastOutputPtr blastOP = NULL;
	int evt = 1;
	IterationPtr itFor = NULL;
	IterationPtr itRev = NULL;
	ShortReadPtr seqFor = NULL;
	ShortReadPtr seqRev = NULL;
	gzFile fp = NULL;
	gzFile fp2 = NULL;
	int inter = 1; // NEED TO PUT IT IN OPTION

	if (argc < 4)
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

	if (parseDict(argv[2]))
		ERROR("Error while reading the index base", EXIT_FAILURE)

	fp = gzopen(argv[3], "r");
	if (fp == NULL)
		ERROR("Unable to open the FastQ\n", EXIT_FAILURE)
	
	if (argc == 5 && !inter)	
	{
		fp2 = gzopen(argv[4], "r");
		if (fp2 == NULL)
			ERROR("Unable to open the second FastQ\n", EXIT_FAILURE)
	}

	while (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "Iteration"))
	{
		itFor = parseIteration(reader);
		seqFor = shortReadNext(fp);
		
		if (fp2 != NULL || inter)
		{
			itRev = parseIteration(reader);
			if (inter)
				seqRev = shortReadNext(fp);
			else
				seqRev = shortReadNext(fp2);
			printSAM(itFor, itRev, seqFor, seqRev);
			deallocIteration(itRev);
			shortReadFree(seqRev);
		}

		else
			printSAM(itFor, NULL, seqFor, NULL);

		deallocIteration(itFor);
		shortReadFree(seqFor);	
	}
	
	gzclose(fp);
	if (fp2 != NULL)
		gzclose(fp2);
		
	deallocBlastOutput(blastOP);
	
	xmlFreeTextReader(reader);
	xmlCleanupCharEncodingHandlers();
	xmlDictCleanup();
	
	return EXIT_SUCCESS;
}





