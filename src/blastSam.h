
#ifndef SAM_H_INCLUDED
#define SAM_H_INCLUDED

typedef struct AppParam
{
	int inter;
	int minLen;
	char* blastOut;
	char* db;
	char* fastq1;
	char* fastq2;
}AppParam, *AppParamPtr;

int blastToSam(AppParamPtr app);

#endif // SAM_H_INCLUDED
