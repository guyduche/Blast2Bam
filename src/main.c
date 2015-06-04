
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "blastSam.h"
#include "utils.h"

int main(int argc, char** argv)
{
	AppParamPtr app = safeCalloc(1, sizeof(AppParam));
	int ret = 0, c = 0, option_index = 0;
	static struct option long_options[] = 
	{
		{"interleaved", no_argument, 0, 'p'},
		{"minAlignLength", required_argument, 0, 'W'},
		{"posOnChr", no_argument, 0, 'z'},
		{"readGroup", required_argument, 0, 'R'},
		{0, 0, 0, 0}
	};
	
	while ((c = getopt_long(argc, argv,"pR:W:z", long_options, &option_index)) != -1)
	{
		switch (c)
		{
			case 'p': app->inter = 1; break;
			case 'R': app->readGroup = optarg; break;
			case 'W': app->minLen = (int) strtol(optarg, NULL, 10); break;
			case 'z': app->posOnChr = 1; break;
			case ':': fprintf(stderr, "Option -%c requires an argument\n", optopt); return EXIT_FAILURE; break;
			case '?': fprintf(stderr, "Option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
		}
	}
	
	if (argc < optind + 2)
		return EXIT_FAILURE;
	
	app->blastOut = argv[optind];
	app->db = argv[optind + 1];
	app->fastq1 = argv[optind + 2];
	
	if (argc == optind + 4 && !app->inter)
		app->fastq2 = argv[optind + 3];
	
	ret = blastToSam(app);
	
	free(app);
	
	return ret;
}
