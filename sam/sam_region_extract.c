#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "reader.h"
#include "region.h"

/*
int main(int argc, char * argv[]) {
	char buffer[512];
	sam_file_t * samfile = sam_open(argv[1], buffer, 512);
	if(!samfile) return 1;
	sam_line_t * samline = NULL;
	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return 1;
	}
	int extension = atoi(argv[2]);
	samline = sam_readline(samfile, samline, buffer, 512, extension);
	region->start = samline->pos;
	region->end = samline->pos + samline->read_len;
	strcpy(region->chrom, samline->chromosome);
	region = get_region(samfile, samline, region, 0, extension, buffer, 512);
	printf("%d\t%d\n", region->start, region->end);
	region = get_region(samfile, samline, region, 1, extension, buffer, 512);
	printf("%d\t%d\n", region->start, region->end);
	region = get_region(samfile, samline, region, 1, extension, buffer, 512);
	printf("%d\t%d\n", region->start, region->end);

	return 0;
}
*/

region_t * get_region(sam_file_t * samfile, sam_line_t * samline, region_t * region, int new_region, int extension, char * buf, int buf_len) {
	if(new_region) {
		//Starting of a new region, store chromosome
		strcpy(region->chrom, samline->chromosome);
		region->start = samline->pos;
		region->end = samline->pos + samline->read_len;
		new_region = 0;
	} else {
		samline = sam_readline(samfile, samline, buf, buf_len, extension);
	}

	if(samline->flag != 4) {
		if(samline->pos > region->end || samline->pos < region->start) {
			//out of region bound, should return start/end and store this pos
			new_region = 1;
			return region;
		}

		//pos is whitin region bounds
		region->end = samline->pos + samline->read_len;
	}

	return get_region(samfile, samline, region, new_region, extension, buf, buf_len);
}
