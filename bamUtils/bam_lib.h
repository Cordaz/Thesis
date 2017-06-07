#ifndef BAM_LIB_H
	#define BAM_LIB_H

	#include <htslib/hts.h>
	#include <htslib/sam.h>
	#include "genome.h"

	#define REG_COMPLETE 0
	#define CONTINUE 1
	#ifndef SEQUENCE_BUFFER
		#define SEQUENCE_BUFFER 1048576
	#endif

	/********** STRUCT **********/
	/*
	 * BAM file structure
	 * @field in pointer to BAM/SAM file (htslib)
	 * @field header pointer to BAM header struct (htslib)
	 * @field aln pointer to BAM alignment struct (htslib)
	 */
	typedef struct myBam_s {
		samFile * in;
		bam_hdr_t * header;
		bam1_t * aln;
	} myBam_t;

	/*
	 * Seqence structure
	 * @field seq nucleotides sequence
	 * @field len length of the sequence
	 */
	typedef struct sequence_s {
		char seq[SEQUENCE_BUFFER];
		int len;
	} sequence_t;


	/*
	 * Opens a BAM/SAM file, reads the header and initialize an alignment
	 * @param bam_path
	 * @return a BAM struct
	 */
	myBam_t * myBam_open(char *);


	/*
	 * Close a BAM file and frees the struct
	 */
	void myBam_close(myBam_t *);

	/*
	 * Requires both mybam and region to be well initialized, skip control
	 * @param myBam BAM struct (see definition above)
	 * #param region region struct (see definition above)
	 * @param extension size of extension (otherwise use original size)
	 * @param *status is 0 if region is completed, 1 if is still under construction, -1 if bam is EOF
	 * @return region or NULL if error occurred
	 *
	 */
	region_t * get_next_region(myBam_t *, region_t *, int, int *);

	/*
	 * Get the region before (i.e. actually before if pos strand, otherwise the after) of the given region.
	 * @param region region struct of which to find the before
	 * @param before_region the initialized region struct to build the to return region
	 * @param chrom_info chromosome_info struct to guarantee that the region would be within the chromosome border
	 * @return before_region
	 *
	 */
	region_t * get_before_region(region_t *, region_t *, chromosomes_info_t *);


	/*
	 * Requires both mybam and region to be well initialized, skip control
	 * @param myBam BAM struct (see definition above)
	 * @param region region struct (see definition above)
	 * @param extension size of extension (otherwise use original size)
	 * @param *status is 0 if region is completed, 1 if is still under construction, -1 if bam is EOF
	 * @return region or NULL if error occurred
	 *
	 */
	region_t * get_next_region_overlap(myBam_t *, region_t *, int, int *);


	/*
	 * Requires genome initialized and region to be allocated
	 * @param genome the genome struct initialized
	 * @region region from which to extract the sequence
	 * @param sequence sequence struct to be used
	 * @return updated sequence struct
	 */
	sequence_t * get_sequence(genome_t *, region_t *, sequence_t *);
#endif
