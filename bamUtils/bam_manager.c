#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <malloc.h>

#include "bam_manager.h"

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
