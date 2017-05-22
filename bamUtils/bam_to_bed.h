#ifndef BAM_TO_BED_H
	#define BAM_TO_BED_H

	#include "bam_to_bed.h"
	#include "bam_to_region.h"
	#include "genome.h"
	#include "bam_manager.h"

	/*
	 * bam_path, bed_path, extension size, boolean region (either consider single read or complete region)
	 *
	 */
	int bam_to_bed(char *, char *, int, int);

#endif
