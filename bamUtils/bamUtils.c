#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../utilities/argparse.h"
#include "bam_lib.h"
#include "genome.h"

#ifndef BUFFER
	#define BUFFER 512
#endif

static const char* const usage[] = {
	"bamUtils [options] (-b|-fg <path/to/genome/dir/>|-bfg <path/to/genome/dir/>) -i <in.bam>|<in.sam>",
	NULL
};

int main(int argc, const char * argv[]) {
	const char * bam_file = NULL;
	const char * genome_dir = NULL;
	int bed_arg = 0;
	int fa_arg = 0;
	int l_arg = 0;
	int region_arg = 0;

	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Required"),
		OPT_STRING('i', "input", &bam_file, "BAM/SAM input file"),
		OPT_GROUP("Select either b, f or both (f requires g arg)"),
		OPT_BOOLEAN('b', "bed", &bed_arg, "generate BED file"),
		OPT_BOOLEAN('f', "fasta", &fa_arg, "generate FASTA file"),
		OPT_STRING('g', "genome", &genome_dir, "path/to/genome/dir/"),
		OPT_GROUP("Optional"),
		OPT_INTEGER('l', "length", &l_arg, "length of the extension to be considered (default no extension)"),
		OPT_BOOLEAN('r', "region", &region_arg, "consider overlapping region as unique region (i.e. multiple reads are considered as one region if partially overlapping)"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nConvert a BAM/SAM file to a region (BED file) or sequence (FASTA file)", "\nCan specify an extension (otherwise the original length is preserved) and if to consider overlapping reads as single region");

	argc = argparse_parse(&argparse, argc, argv);

	if(bam_file == NULL) {
		fprintf(stdout, "[ERROR] must specify a BAM/SAM input file\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(!bed_arg && !fa_arg) {
		fprintf(stdout, "[ERROR] must specify an output type\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(fa_arg && genome_dir == NULL) {
		fprintf(stdout, "[ERROR] must specify the genome directory for creating the FASTA file\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	char support[BUFFER+1];
	char support2[BUFFER+1];
	char * extension;
	char bed_file[BUFFER+1];
	strncpy(bed_file, bam_file, BUFFER);
	extension = strrchr(bed_file, '.');
	strcpy(support, "");
	strcpy(support2, "");
	if(l_arg > 0) {
		snprintf(support2, BUFFER, "_ext%d", l_arg);
	}
	strcat(support, support2);
	strcpy(support2, "");
	if(region_arg) strcpy(support2, "_region");
	strcat(support, support2);
	//printf("%s\n", support);
	strcpy(extension, support);
	strcat(extension, ".bed");
	char fa_file[BUFFER+1];
	strncpy(fa_file, bam_file, BUFFER);
	extension = strrchr(fa_file, '.');
	strcpy(extension, support);
	strcat(extension, ".fa");


	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return -1;
	}

	myBam_t * myBam = myBam_open((char*)bam_file);
	chromosomes_info_t * chrom_info;
	genome_t * genome = NULL;
	sequence_t * sequence = NULL;

	if(fa_arg) {
		if(!(chrom_info = chromosome_info_init(myBam->header->target_name, myBam->header->target_len, myBam->header->n_targets))) {
			fprintf(stdout, "[ERROR] failed to load chromosomes info\n");
			return 1;
		}
		genome = genome_init((char*)genome_dir, chrom_info);
	}


	FILE * fa_fp;
	if(fa_arg) {
		if(!(fa_fp = fopen(fa_file, "w+"))) {
			fprintf(stdout, "[ERROR] can't open %s\n", fa_file);
		}
	}

	FILE * bed_fp;
	if(bed_arg) {
		if(!(bed_fp = fopen(bed_file, "w+"))) {
			fprintf(stdout, "[ERROR] can't open %s\n", bed_file);
		}
		fprintf(bed_fp, "chrom\tchromStart\tchromEnd\n");
	}

	int status = REG_COMPLETE;

	if(region_arg) { //TO_REGION
		status = sam_read1(myBam->in, myBam->header, myBam->aln);
		if(status <= 0) {
			fprintf(stdout, "[ERROR] unexpected EOF\n");
			return -2;
		}
		status = REG_COMPLETE;

		region = get_next_region_overlap(myBam, region, l_arg, &status);
		while( status != EOF) {
			if(status == REG_COMPLETE) {
				if(fa_arg) {
					if(!(sequence = get_sequence(genome, region, sequence))) {
						return 1;
					}
					fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
				}
				if(bed_arg) fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end-1);
			}
			region = get_next_region_overlap(myBam, region, l_arg, &status);
			//printf("%d\n", status);
		}
		if(fa_arg) {
			if(!(sequence = get_sequence(genome, region, sequence))) {
				return 1;
			}
			fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
		}
		if(bed_arg) fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end-1);
	} else {
		region = get_next_region(myBam, region, l_arg, &status);
		while( status != EOF ) {
			//printf("%s:%d-%d\n", region->chromosome, region->start, region->end);
			if(status == REG_COMPLETE) {
				if(fa_arg) {
					if(!(sequence = get_sequence(genome, region, sequence))) {
						return 1;
					}
					fprintf(fa_fp, ">%s:%d\n%s\n", region->chromosome, region->start, sequence->seq);
				}
				if(bed_arg) fprintf(bed_fp, "%s\t%d\t%d\n", region->chromosome, region->start, region->end-1);
			}
			region = get_next_region(myBam, region, l_arg, &status);
		}
	}

	if(fa_arg) fclose(fa_fp);
	if(bed_arg) fclose(bed_fp);

	myBam_close(myBam);

	return 0;

}
