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


int main(int argc, char * argv[]) {
	genome_t * genome = NULL;
	sequence_t * sequence = NULL;
	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return -1;
	}

	//TEST CASE - OPEN A BAM - NOT NEEDED AFTER
	myBam_t * myBam = myBam_open(argv[1]);
	chromosomes_info_t * chrom_info;
	if(!(chrom_info = chromosome_info_init(myBam->header->target_name, myBam->header->target_len, myBam->header->n_targets))) {
		fprintf(stdout, "[ERROR] failed to load chromosomes info\n");
		return 1;
	}
	genome = genome_init(argv[2], chrom_info);

	FILE * fa_fp;
	if(!(fa_fp = fopen(argv[3], "w+"))) {
		fprintf(stdout, "[ERROR] can't open %s\n", argv[3]);
	}


	if(atoi(argv[5])) { //TO_REGION
		int status;
		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) {
			fprintf(stdout, "[ERROR] unexpected EOF\n");
			return -2;
		}

		while(region = get_next_region_overlap(myBam, region, atoi(argv[4]), 1)) {
			//printf("Region: %s:%d-%d\n", region->chromosome, region->start, region->end);
			if(!(sequence = get_sequence(genome, region, sequence))) {
				return 1;
			}
			fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start+1, sequence->seq);
		}
	} else {
		while(region = get_next_region(myBam, region, atoi(argv[4]))) {
			//printf("%s:%d-%d\n", region->chromosome, region->start, region->end);

			if(!(sequence = get_sequence(genome, region, sequence))) {
				return 1;
			}

			fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start+1, sequence->seq);

		}
	}

	fclose(fa_fp);

	return 0;
}
