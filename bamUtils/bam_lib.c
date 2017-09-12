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
			if(end > myBam->header->target_len[myBam->aln->core.tid]) end = myBam->header->target_len[myBam->aln->core.tid]; //Out of chromosome
		}
		region->strand = '+';
	} else if(myBam->aln->core.flag == 16) {
		if(extension) {
			start = end - extension;
			if(start < 0) start = 0; //Out of chromosome
		}
		region->strand = '-';
	}

	if( strcmp(region->chromosome, myBam->header->target_name[myBam->aln->core.tid]) == 0 ) { //Remove duplicates
		if(start <= region->start) {
			*status = CONTINUE;
			return region;
		}
	}

	strncpy(region->chromosome, myBam->header->target_name[myBam->aln->core.tid], 6);
	region->chrom_index = myBam->aln->core.tid;
	region->start = start;
	if(end > myBam->header->target_len[myBam->aln->core.tid]) {
		end = myBam->header->target_len[myBam->aln->core.tid];
	}

	region->end = end;

	*status = REG_COMPLETE;
	return region;
}

region_t * get_before_region(region_t * region, region_t * before_region, chromosomes_info_t * chrom_info, int skip) {
	int region_size = region->end - region->start;
	if(region->strand == '+') {
		before_region->start = region->start - region_size - skip;
		if(before_region->start < 0) before_region->start = 1;
		before_region->end = region->start - skip;
	} else {
		before_region->start = region->end + skip;
		before_region->end = region->end + region_size + skip;
		if(before_region->end > chrom_info->sizes[region->chrom_index]) before_region->end = chrom_info->sizes[region->chrom_index];
	}

	strcpy(before_region->chromosome, region->chromosome);
	before_region->chrom_index = region->chrom_index;

	return before_region;
}

region_t * get_after_region(region_t * region, region_t * after_region, chromosomes_info_t * chrom_info, int skip) {
	int region_size = region->end - region->start;
	if(region->strand == '+') {
		after_region->start = region->end + skip;
		if(after_region->start < 0) after_region->start = 1;
		after_region->end = after_region->start + region_size;
	} else {
		after_region->end = region->start - skip;
		after_region->start = region->start - region_size - skip;
	}
	if(after_region->end > chrom_info->sizes[region->chrom_index]) after_region->end = chrom_info->sizes[region->chrom_index];

	strcpy(after_region->chromosome, region->chromosome);
	after_region->chrom_index = region->chrom_index;

	return after_region;
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
		region->chrom_index = myBam->aln->core.tid;
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
		if(!(genome = genome_load_chromosome(genome, region->chromosome, region->chrom_index)) ) {
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
