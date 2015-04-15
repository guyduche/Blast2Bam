#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parseXML.h"
#include "cigar.h"
#include "macro.h"

#define CSTRMACRO(TAG, FUNCT) { \
		COUNTCHARMACRO(FUNCT); \
		cigarStr[*sizeCStr-2] = count; \
		cigarStr[*sizeCStr-1] = TAG;}
		
#define COUNTCHARMACRO(FUNCT) { \
		count = 1; \
		pos++; \
		while (pos < (hsp->hsp_align_len) && FUNCT) \
		{ \
			count++; \
			pos++; \
		} \
		pos--;}

int *cigarStrBuilding(int *cigarStr, Hsp *hsp, int queryLength, int *sizeCStr)
{
	int pos = 0;
	int count = 0;
	*sizeCStr = 0;
	
	if (queryLength > (hsp->hsp_align_len) && (hsp->hsp_query_from) > 1)
	{
		*sizeCStr += 2;
		TABBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		cigarStr[0] = hsp->hsp_query_from - 1;
		cigarStr[1] = 'S';
	}
	
	for (pos = 0; pos < (hsp->hsp_align_len); pos++)
	{	
		*sizeCStr += 2;
		TABBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		
		if (hsp->hsp_hseq[pos] == '-')
			CSTRMACRO('I', (hsp->hsp_hseq[pos] == '-'))
		
		else if (hsp->hsp_qseq[pos] == '-')
			CSTRMACRO('D', (hsp->hsp_qseq[pos] == '-'))
			
		else if (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos])
			CSTRMACRO('=', (hsp->hsp_hseq[pos] == hsp->hsp_qseq[pos]))
		
		else
			CSTRMACRO('X', (hsp->hsp_hseq[pos] != '-' && hsp->hsp_qseq[pos] != '-' && hsp->hsp_hseq[pos] != hsp->hsp_qseq[pos]))
		
		if (cigarStr[*sizeCStr-2] >= 100 && cigarStr[*sizeCStr-1] == 'D')
			cigarStr[*sizeCStr-1] = 'N';
		
	}
	
	if (queryLength > (hsp->hsp_align_len) && (queryLength - hsp->hsp_query_to) > 1)
	{
		*sizeCStr += 2;
		TABBUILDINGMACRO(cigarStr, *sizeCStr, int, NULL)
		cigarStr[*sizeCStr-2] = queryLength - hsp->hsp_query_to;
		cigarStr[*sizeCStr-1] = 'S';
	}
	
	return cigarStr;
}







