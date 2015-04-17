
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "macro.h"

int main(int argc, char **argv)
{
	FILE *reader;
	char *str = NULL;
	int c;
	int countStr = 0;
	int countSpace = 0;
	
	if (argc != 2)
		ERROR("Wrong number of arguments\n", EXIT_FAILURE)

	reader = fopen(argv[1], "r");

	if (reader == NULL)
		ERROR("Error while opening the file\n", EXIT_FAILURE)

	do
	{
		if (countSpace > 2 || !countSpace)
			do
			{
				c = fgetc(reader);
			} while (c != '\n' && c != EOF);

		countStr = 0;
		countSpace = 0;
		str = NULL;
		c = fgetc(reader);
		countStr++;

		while (c != '\n' && c != EOF && countSpace <= 2)
		{
			if (c == '\t')
				countSpace++;			
			
			TABBUILDINGMACRO(str, countStr+1, char, EXIT_FAILURE)
			str[countStr-1] = (char) c;
			c = fgetc(reader);
			countStr++;
		}

		TABBUILDINGMACRO(str, countStr, char, EXIT_FAILURE)
		str[countStr-1] = '\0';

		fprintf(stdout, "%s", str);
		free(str);
		if (c != EOF)
			fprintf(stdout, "\n");

	} while (c != EOF);
	
	fclose(reader);

	return EXIT_SUCCESS;
}



