#ifndef DEBUG_H
#define DEBUG_H
#include <stdio.h>

#define DEBUG(FORMAT,ARGS...) \
	do {\
	fprintf(stderr,"[%s:%d]:",__FILE__,__LINE__); \
	fprintf(stderr,FORMAT,ARGS);\
	fputc('\n',stderr);\
	}while(0)


#endif

