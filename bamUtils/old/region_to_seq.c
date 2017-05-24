#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "region_to_seq.h"

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

/*
int bam_to_fa(char * bam_file, char * fa_file, char * genome_dir, int extension, int to_region) {
	genome_t * genome = NULL;
	sequence_t * sequence = NULL;
	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return -1;
	}

	myBam_t * myBam = myBam_open(bam_file);
	chromosomes_info_t * chrom_info;
	if(!(chrom_info = chromosome_info_init(myBam->header->target_name, myBam->header->target_len, myBam->header->n_targets))) {
		fprintf(stdout, "[ERROR] failed to load chromosomes info\n");
		return 1;
	}
	genome = genome_init(genome_dir, chrom_info);

	FILE * fa_fp;
	if(!(fa_fp = fopen(fa_file, "w+"))) {
		fprintf(stdout, "[ERROR] can't open %s\n", argv[3]);
	}

	int status = REG_COMPLETE;

	if(to_region) { //TO_REGION
		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) {
			fprintf(stdout, "[ERROR] unexpected EOF\n");
			return -2;
		}
		status = REG_COMPLETE;

		region = get_next_region_overlap(myBam, region, extension, &status);
		while( status != EOF) {
			if(status == REG_COMPLETE) {
				if(!(sequence = get_sequence(genome, region, sequence))) {
					return 1;
				}
				//printf("%s:%d-%d\n", region->chromosome, region->start, region->end);
				fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
			}
			region = get_next_region_overlap(myBam, region, extension, &status);
			//printf("%d\n", status);
		}
		if(!(sequence = get_sequence(genome, region, sequence))) {
			return 1;
		}
		fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
	} else {
		region = get_next_region(myBam, region, extension, &status);
		while( status != EOF ) {
			//printf("%s:%d-%d\n", region->chromosome, region->start, region->end);
			if(status == REG_COMPLETE) {
				if(!(sequence = get_sequence(genome, region, sequence))) {
					return 1;
				}
				fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
			}
			region = get_next_region(myBam, region, extension, &status);
		}
	}

	fclose(fa_fp);

	myBam_close(myBam);

	return 0;
}
*/
