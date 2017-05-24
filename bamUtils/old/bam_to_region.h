#ifndef BAM_TO_REGION_H
	#define BAM_TO_REGION_H

	#include <htslib/hts.h>
	#include <htslib/sam.h>
	#include "genome.h"
	#include "bam_manager.h"

	#define REG_COMPLETE 0
	#define CONTINUE 1

	/*
	 * Requires both mybam and region to be well initialized, skip control
	 *
	 */
	region_t * get_next_region(myBam_t *, region_t *, int, int *);

	region_t * get_next_region_overlap(myBam_t *, region_t *, int, int *);

#endif
