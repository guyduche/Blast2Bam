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
#include <getopt.h>
#include <string.h>
#include "blastSam.h"
#include "utils.h"

static void usage(FILE* out)
{
    fprintf(out, "Blast2Bam. Last compilation: %s at %s.\n\n", __DATE__, __TIME__);
    fprintf(out, "Usage: blast2bam [options] <Blast XML output> <reference sequence dictionary> <FastQ_1> [FastQ_2]\n\n");
    fprintf(out, "Options:\n");
    fprintf(out, " --output         | -o FILE       Output file (default: stdout)\n");
    fprintf(out, " --interleaved    | -p            Interleaved data\n");
    fprintf(out, " --readGroup      | -R STR        Read group header line '@RG\\tID:foo'\n");
    fprintf(out, " --minAlignLength | -W INT        Discard alignments shorter than [INT]\n");
    fprintf(out, " --posOnChr       | -z            Adjust the alignment position to the first position of the reference\n");
    fprintf(out, " --help           | -h            Get help (this screen)\n");
    fprintf(out, "\n\n");
}

int main(int argc, char** argv)
{
    AppParamPtr app = NULL;
    int i, ret = 0, c = 0, option_index = 0;
    FILE* out = NULL;
    static struct option long_options[] =
    {
        {"interleaved", no_argument, 0, 'p'},
        {"minAlignLength", required_argument, 0, 'W'},
        {"posOnChr", no_argument, 0, 'z'},
        {"output", required_argument, 0, 'o'},
        {"readGroup", required_argument, 0, 'R'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    if (argc < 4)
    {
        usage(stderr);
        return EXIT_FAILURE;
    }

    // Create new AppParam
    app = safeCalloc(1, sizeof(AppParam));

    // Get program options
    while ((c = getopt_long(argc, argv,"ho:pR:W:z", long_options, &option_index)) != -1)
    {
        switch (c)
        {
            case 'o': app->out = out = safeFOpen(optarg, "a"); break;                   // Program output
            case 'p': app->inter = 1; break;                                            // Interleaved
            case 'R': app->readGroup = optarg; break;                                   // Read group
            case 'W': app->minLen = (int) strtol(optarg, NULL, 10); break;              // Minimum alignment length
            case 'z': app->posOnChr = 1; break;                                         // Alignment position adjusted to the position of the reference
            case 'h': usage(stderr); return EXIT_SUCCESS; break;                        // Help
            case ':': fprintf(stderr, "Option -%c requires an argument\n", optopt); return EXIT_FAILURE; break;
            case '?': fprintf(stderr, "Option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
        }
    }

    // Set default output to stdout
    if (out == NULL) app->out = stdout;

    if (argc < optind + 3)
    {
        usage(stderr);
        return EXIT_FAILURE;
    }

    // PG line of SAM header
    app->pg_line = safeStrdup("@PG\tID:Blast2Bam\tPN:Blast2Bam\tVN:0.1\tCL:");
    for (i = 0; i < argc; i++)
    {
        if (i > 0) safeStrAppend(app->pg_line, " ");
        safeStrAppend(app->pg_line, argv[i]);
    }

    app->blastOut = argv[optind];       // Get Blast output
    app->db = argv[optind + 1];         // Get the reference sequence dictionary
    app->fastq1 = argv[optind + 2];     // Get the first FastQ

    // A second FastQ file is present
    if (argc == optind + 4 && !app->inter)
    {
        app->fastq2 = argv[optind + 3]; // Get the second FastQ
        if (!strcmp(app->fastq1, app->fastq2))
            ERROR("FastQ_1 and FastQ_2 must be different\n", EXIT_FAILURE);
    }

    // Call to the main function of Blast2Sam
    ret = blastToSam(app);

    if (out != NULL) fclose(out);
    if (app->readGroupID != NULL)
        free(app->readGroupID);
    free(app->pg_line);
    free(app);
    return ret;
}
