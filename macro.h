#ifndef MACRO_H_INCLUDED
#define MACRO_H_INCLUDED

#define ERROR(MESSAGE, FUNCT, ENDING) { \
		fprintf(stderr, MESSAGE); \
		FUNCT; \
		return ENDING;}

#define TABBUILDINGMACRO(TABNAME, NBALLOC, TYPE, ENDING) TABNAME = (TYPE*) realloc(TABNAME, (NBALLOC * sizeof(TYPE))); \
		if (TABNAME == NULL) \
			ERROR ("Bad memory allocation of TABNAME", NULL, ENDING)

#endif // MACRO_H_INCLUDED
