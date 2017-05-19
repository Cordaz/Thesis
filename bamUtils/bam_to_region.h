#ifndef BAM_TO_REGION_H
	#define BAM_TO_REGION_H

	#include <htslib/hts.h>
	#include <htslib/sam.h>

	typedef struct myBam_s {
		samFile * in;
		bam_hdr_t * header;
		bam1_t * aln;
	} myBam_t;

	typedef struct region_s {
		char chromosome[6];
		int start;
		int end;
	} region_t;

	myBam_t * myBam_open(char *);

	void myBam_close(myBam_t *);

	/*
	 * Requires both mybam and region to be well initialized, skipp control
	 *
	 */
	region_t * get_next_region(myBam_t *, region_t *, int);

	region_t * get_next_region_overlap(myBam_t *, region_t *, int, int);

#endif
