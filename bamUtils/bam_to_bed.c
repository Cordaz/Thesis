#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

//#include "bam_to_bed.h"
#include "bam_to_region.h"
#include "genome.h"

int bam_to_bed(char * bam, char * bed, int extension, int to_region) {

	myBam_t * myBam = myBam_open(bam);

	FILE * bed_fp;
	if(!(bed_fp = fopen(bed, "w+"))) {
		fprintf(stdout, "[ERROR] can't open %s\n", bed);
		return 1;
	}

	fprintf(bed_fp, "chrom\tchromStart\tchromEnd\n");

	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return -1;
	}
	int status = REG_COMPLETE;

	if(to_region) {
		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) {
			fprintf(stdout, "[ERROR] unexpected EOF\n");
			return -2;
		}
		status = REG_COMPLETE;

		region = get_next_region_overlap(myBam, region, extension, &status);
		while( status != EOF ) {
			if(status == REG_COMPLETE) {
				//printf("%s\t%d\t%d\n", region->chromosome, region->start, region->end);
				fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end);
			}
			//printf("%d\n", status);
			region = get_next_region_overlap(myBam, region, extension, &status);
		}
		fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end);
	} else {
		region = get_next_region(myBam, region, extension, &status);
		while( status != EOF) {
			if(status == REG_COMPLETE) {
				fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end);
			}
			region = get_next_region(myBam, region, extension, &status);
		}
	}


	fclose(bed_fp);

	myBam_close(myBam);

	return 0;
}


int main(int argc, char * argv[]) {
	return bam_to_bed(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
}
