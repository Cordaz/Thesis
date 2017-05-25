#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "bam_lib.h"

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

	if(start < 0) start=0;

	if( strcmp(region->chromosome, myBam->header->target_name[myBam->aln->core.tid]) == 0 ) { //Remove duplicates
		if(start <= region->start) {
			*status = CONTINUE;
			return region;
		}
	}

	strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 6);
	region->start = start;
	if(end > myBam->header->target_len[myBam->aln->core.tid]) {
		end = myBam->header->target_len[myBam->aln->core.tid];
	}

	region->end = end;

	*status = REG_COMPLETE;
	return region;
}

region_t * get_next_region_overlap(myBam_t * myBam, region_t * region, int extension, int * status) {
	int file_status;

	if(*status != REG_COMPLETE) {
		file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(file_status <= 0) {
			*status = EOF;
			return region;
		}
	}

	if(myBam->aln->core.flag != 0 && myBam->aln->core.flag != 16) {
		if(*status == REG_COMPLETE) {
			file_status = sam_read1(myBam->in, myBam->header, myBam->aln);
			if(file_status <= 0) {
				*status = EOF;
				return region;
			}
			return region;
		}
		*status = CONTINUE;
		return region;
	}

	int startpos, endpos;
	startpos = myBam->aln->core.pos; //1 based
	endpos = bam_endpos(myBam->aln);

	if(myBam->aln->core.flag == 0) {
		if(extension) {
			endpos = startpos + extension;
		}
	} else if(myBam->aln->core.flag == 16) {
		if(extension) {
			startpos = endpos - extension;
			if(startpos < 0) startpos = 0;
		}
	}

	if(endpos > myBam->header->target_len[myBam->aln->core.tid]) {
		endpos = myBam->header->target_len[myBam->aln->core.tid];
	}

	if(*status == REG_COMPLETE) {
		strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 6);
		region->start = startpos;
		region->end = endpos;

		*status = CONTINUE;
		return region;
	}

	//A region already exists, check if can be extended.
	if(startpos <= region->end && startpos >= region->start) {
		region->end = endpos;

		*status = CONTINUE;
		return region;
	}
	if(endpos <= region->end && endpos >= region->start) {
		region->start = startpos;

		*status = CONTINUE;
		return region;
	}

	//Can't be extended, return the region created up to now
	*status = REG_COMPLETE;
	return region;
}

sequence_t * get_sequence(genome_t * genome, region_t * region, sequence_t * sequence) {
	if(!region) return NULL;

	if( !genome->chromosome || strcmp(region->chromosome, genome->chromosome->name) != 0 ) {
		if(!(genome = genome_load_chromosome(genome, region->chromosome)) ) {
			fprintf(stdout, "[ERROR] can't load correctly %s\n", region->chromosome);
			return NULL;
		}
	}

	if(!sequence) {
		if(!(sequence = (sequence_t*)malloc(sizeof(sequence_t)))) {
			fprintf(stdout, "[ERROR] can't allocate\n");
			return NULL;
		}
	}

	int i, j;
	for(i=region->start, j=0; i<region->end && j < SEQUENCE_BUFFER; i++, j++) {
		sequence->seq[j] = genome->chromosome->seq[i];
	}
	sequence->seq[j] = '\0';
	sequence->len = j;

	return sequence;
}
