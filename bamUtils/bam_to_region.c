#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "bam_to_region.h"

region_t * get_next_region(myBam_t * myBam, region_t * region, int extension, int * status) {
	int file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
	if(file_status <= 0) {
		*status = EOF;
		return region;
	}

	if(myBam->aln->core.flag != 0 && myBam->aln->core.flag != 16) {
		*status = CONTINUE;
		return region;
	}

	int start = 0;
	int end = 0;

	start = myBam->aln->core.pos;
	end = bam_endpos(myBam->aln);


	if(myBam->aln->core.flag == 0) {
		if(extension) {
			end = start + extension;
		}
	} else if(myBam->aln->core.flag == 16) {
		if(extension) {
			start = end - extension;
		}
	}


	if( strcmp(region->chromosome, myBam->header->target_name[myBam->aln->core.tid]) == 0 ) { //Remove duplicates
		if(start <= region->start) {
			*status = CONTINUE;
			return region;
		}
	}

	strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 6);
	region->start = start;
	if(end > myBam->header->target_len[myBam->aln->core.tid]) {
		*status = 1;
		return region;
	}

	region->end = end;

	*status = REG_COMPLETE;
	return region;
}

region_t * get_next_region_overlap(myBam_t * myBam, region_t * region, int extension, int * new_region, int * status) {
	int file_status;

	if(myBam->aln->core.flag != 0 && myBam->aln->core.flag != 16) {
		//if(*new_region) {
			file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
			if(file_status <= 0) {
				*status = EOF;
				return region;
			}
		//}
		*status = CONTINUE;
		return region;
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


	if(endpos > myBam->header->target_len[myBam->aln->core.tid]) {
		endpos = myBam->header->target_len[myBam->aln->core.tid];
	}

	if(*new_region) {
		strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 6);
		region->start = startpos;
		region->end = endpos;

		file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(file_status <= 0) {
			*status = EOF;
			return region;
		}

		*status = CONTINUE;
		*new_region = 0;
		return region;
	}

	//A region already exists, check if can be extended.
	if(startpos <= region->end && startpos >= region->start) {
		region->end = endpos;

		file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(file_status <= 0) {
			*status = EOF;
			return region;
		}

		*status = CONTINUE;
		*new_region = 0;
		return region;
	}

	//Can't be extended, return the region created up to now
	*status = REG_COMPLETE;
	*new_region = 1;
	return region;
}
