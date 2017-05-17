#ifndef SAM_REGION_H
	#define SAM_REGION_H

	#include "reader.h"

	typedef struct region_s {
		int start, end;
		char chrom[3];
	} region_t;

	region_t * get_region(sam_file_t *, sam_line_t *, region_t *, int, int, char *, int);

#endif
