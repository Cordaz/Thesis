#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "bam_to_region.h"


myBam_t * myBam_open(char * bamfile) {
	myBam_t * myBam;
	if( !(myBam = (myBam_t*)malloc(sizeof(myBam_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate \n");
		return NULL;
	}

	myBam->in = sam_open(bamfile, "r");
	myBam->header = sam_hdr_read(myBam->in);
	myBam->aln = bam_init1();

	return myBam;
}

void myBam_close(myBam_t * myBam) {
	bam_destroy1(myBam->aln);
	bam_hdr_destroy(myBam->header);
	sam_close(myBam->in);
	free(myBam);
}

region_t * get_next_region(myBam_t * myBam, region_t * region, int extension) {
	int status = sam_read1(myBam->in, myBam->header, myBam->aln);
	if(status <= 0) return NULL;

	if(myBam->aln->core.flag == 4) {
		return get_next_region(myBam, region, extension);
	}

	strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 5);
	region->start = myBam->aln->core.pos; //0 based
	region->end = bam_endpos(myBam->aln);

	if(myBam->aln->core.flag == 0) {
		if(extension) {
			region->end = region->start + extension;
		}
	} else if(myBam->aln->core.flag == 16) {
		if(extension) {
			region->start = region->end - extension;
		}
	}

	return region;
}

region_t * get_next_region_overlap(myBam_t * myBam, region_t * region, int extension, int new_region) {
	int status;

	if(myBam->aln->core.flag == 4) {
		if(new_region) {
			status = sam_read1(myBam->in, myBam->header, myBam->aln);
			if(status <= 0) return NULL;
		}
		return get_next_region_overlap(myBam, region, extension, new_region);
	}

	int startpos, endpos;
	startpos = myBam->aln->core.pos; //0 based
	endpos = bam_endpos(myBam->aln);

	if(myBam->aln->core.flag == 0) {
		if(extension) {
			endpos = startpos + extension;
		}
	} else if(myBam->aln->core.flag == 16) {
		if(extension) {
			startpos = endpos - extension;
		}
	}

	if(new_region) {
		strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 5);
		region->start = startpos;
		region->end = endpos;

		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) return NULL;

		return get_next_region_overlap(myBam, region, extension, 0);
	}

	//A region already exists, check if can be extended.
	if(startpos <= region->end && startpos >= region->start) {
		region->end = endpos;

		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) return region;

		return get_next_region_overlap(myBam, region, extension, 0);
	}

	//Can't be extended, return the region created up to now
	return region;
}
