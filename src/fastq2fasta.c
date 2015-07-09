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
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include "utils.h"
#include "shortRead.h"

static void usage(FILE* out)
{
    fprintf(out, "Fastq2fasta. Last compilation: %s at %s.\n\n", __DATE__, __TIME__);
    fprintf(out, "Usage: fastq2fasta [options] <FastQ_1> [FastQ_2]\n\n");
    fprintf(out, "Options:\n");
    fprintf(out, " --output         | -o FILE       Output file (default: stdout)\n");
    fprintf(out, " --help           | -h            Get help (this screen)\n");
    fprintf(out, "\n\n");
}

int main(int argc, char** argv)
{
    gzFile fp = NULL, fp2 = NULL;
    ShortReadPtr reads[2];
    FILE* out = NULL;
    int c = 0, option_index = 0;

    static struct option long_options[] =
    {
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    if (argc < 2)
    {
        usage(stderr);
        return EXIT_FAILURE;
    }

    // Get program options
    while ((c = getopt_long(argc, argv,"ho:", long_options, &option_index)) != -1)
    {
        switch (c)
        {
            case 'o': out = safeFOpen(optarg, "a"); break;                  // Program output
            case 'h': usage(stderr); return EXIT_SUCCESS; break;            // Help
            case ':': fprintf(stderr, "Option -%c requires an argument\n", optopt); return EXIT_FAILURE; break;
            case '?': fprintf(stderr, "Option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
        }
    }

    // Set default output to stdout
    if (out == NULL) out = stdout;

    if (argc < optind)
    {
        usage(stderr);
        return EXIT_FAILURE;
    }

    fp = safeGzOpen(argv[optind], "r");             // Get the first FastQ

    if (argc == optind + 2)
        fp2 = safeGzOpen(argv[optind + 1], "r");    // Get the second FastQ

    reads[0] = shortReadNext(fp);
    while (reads[0] != NULL)
    {
        if (fp2 != NULL)
            reads[1] = shortReadNext(fp2);

        fprintf(out, ">%s\n%s\n", reads[0]->name, reads[0]->seq);
        if (fp2 != NULL)
        {
            fprintf(out, ">%s\n%s\n", reads[1]->name, reads[1]->seq);
            shortReadFree(reads[1]);
        }
        shortReadFree(reads[0]);
        reads[0] = shortReadNext(fp);
    }

    gzclose(fp);
    if (fp2 != NULL) gzclose(fp2);
    if (out != stdout) fclose(out);

    return 0;
}
