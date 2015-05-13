
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "utils.h"
#include "shortRead.h"

int main(int argc, char** argv)
{
	gzFile fp = NULL;
	gzFile fp2 = NULL;
	ShortReadPtr reads[2];
	int inter = 0; // TODO: Put it in option -> get_longopt/getopt (look at Jennifer's code)

	fp = safeGzOpen(argv[1], "r");

	if (argc == 3 && !inter)
		fp2 = safeGzOpen(argv[2], "r");

	reads[0] = shortReadNext(fp);
	while (reads[0] != NULL)
	{
		if (fp2 != NULL)
		{
			if (inter)
				reads[1] = shortReadNext(fp);
			else
				reads[1] = shortReadNext(fp2);
		}

		fprintf(stdout, ">%s\n%s\n", reads[0]->name, reads[0]->seq);
		if (fp2 != NULL)
		{
			fprintf(stdout, ">%s\n%s\n", reads[1]->name, reads[1]->seq);
			shortReadFree(reads[1]);
		}
		shortReadFree(reads[0]);
		reads[0] = shortReadNext(fp);
	}

	gzclose(fp);
	if (fp2 != NULL)
		gzclose(fp2);

	return 0;
}
