
#ifndef MACRO_H_INCLUDED
#define MACRO_H_INCLUDED

#define ERROR(MESSAGE, ENDING) { \
		fprintf(stderr, MESSAGE); \
		return ENDING;}

#define STRBUILDINGMACRO(STRNAME, NBALLOC, TYPE, ENDING) STRNAME = (TYPE*) realloc(STRNAME, (NBALLOC * sizeof(TYPE))); \
		if (STRNAME == NULL) \
			ERROR ("Failed memory allocation of STRNAME", ENDING)

/*
void* _safeMalloc(const char* file,int line,size_t n)
void* _safeCalloc(const char* file,int line,size_t m,size_t n);
void* _safeRealloc(const char* file,int line,void* ptr,size_t n);
FILE* _safeFOPen(...)

#define safeMalloc(n) _safeMalloc(__FILE__,__LINE__,n)
#define safeRealloc(p,n) _safeRealloc(__FILE__,__LINE__,p,n)

*/

#endif // MACRO_H_INCLUDED
