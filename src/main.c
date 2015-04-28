
#include <stdio.h>
#include <stdlib.h>
#include "blastSam.h"


int main(int argc, char** argv)
{
	int ret = 0;
	ret = blastToSam(argc, argv);
	return ret;
}
