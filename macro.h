#ifndef MACRO_H_INCLUDED
#define MACRO_H_INCLUDED

#define ERROR(MESSAGE, ENDING) { \
		fprintf(stderr, MESSAGE); \
		return ENDING;}

#define TABBUILDINGMACRO(TABNAME, NBALLOC, TYPE, ENDING) TABNAME = (TYPE*) realloc(TABNAME, (NBALLOC * sizeof(TYPE))); \
		if (TABNAME == NULL) \
			ERROR ("Bad memory allocation of TABNAME", ENDING)

#endif // MACRO_H_INCLUDED
