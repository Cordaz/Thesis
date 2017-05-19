#ifndef BAM_MNG_H
	#define BAM_MNG_H

	#include "genome.h"

	typedef struct myBam_s {
		samFile * in;
		bam_hdr_t * header;
		bam1_t * aln;
	} myBam_t;

	myBam_t * myBam_open(char *);

	void myBam_close(myBam_t *);

#endif
