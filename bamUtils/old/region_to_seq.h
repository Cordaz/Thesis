#ifndef REGION_TO_SEQ_H
	#define REGION_TO_SEQ_H

	#include "bam_to_region.h"
	#include "genome.h"

	#define SEQUENCE_BUFFER 1048576

	typedef struct sequence_s {
		char seq[SEQUENCE_BUFFER];
		int len;
	} sequence_t;

	/*
	 * Assume genome is fully initialized
	 *
	 */
	sequence_t * get_sequence(genome_t *, region_t *, sequence_t *);

#endif
