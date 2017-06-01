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
	"bamUtils -g <path/to/genome/dir/>) -i <in.bed>|<in.narrowPeak>",
	NULL
};

chromosomes_info_t * open_chromsize(const char *);

int main(int argc, const char * argv[]) {
	const char * bed_file = NULL;
	const char * genome_dir = NULL;

	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Required"),
		OPT_STRING('i', "input", &bed_file, "BED input file"),
		OPT_STRING('g', "genome", &genome_dir, "path/to/genome/dir/"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nConvert a BED-type file to sequence (FASTA file)", "\nIt expects to find genome as chr*.fa and a chromsize file\n");

	argc = argparse_parse(&argparse, argc, argv);

	if(bed_file == NULL) {
		fprintf(stdout, "[ERROR] must specify a BED input file\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(genome_dir == NULL) {
		fprintf(stdout, "[ERROR] must specify the genome directory\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	chromosomes_info_t * chrom_info = open_chromsize(genome_dir);
	if(!chrom_info) return -1;

	genome_t * genome = genome_init((char*)genome_dir, chrom_info);
	if(!genome) return -1;

	FILE * bed_fp;
	if( !(bed_fp = fopen(bed_file, "r")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", bed_file);
		return -1;
	}
	FILE * fa_fp;
	char fa_file[BUFFER];
	strncpy(fa_file, bed_file, BUFFER);
	char * extension = strrchr(fa_file, '.');
	strcpy(extension, ".fa");
	if( !(fa_fp = fopen(fa_file, "w+")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", fa_file);
		return -1;
	}

	char buf[BUFFER];
	char * token;
	region_t * region;
	if( !(region = (region_t*)malloc(sizeof(region_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return -1;
	}
	sequence_t * sequence = NULL;
	//NO HEADER
	fgets(buf, BUFFER, bed_fp);
	while(!feof(bed_fp)) {
		token = strtok(buf, "\t");
		strncpy(region->chromosome, token, 5);
		token = strtok(NULL, "\t");
		region->start = atoi(token);
		token = strtok(NULL, "\t");
		region->end = atoi(token);

		sequence = get_sequence(genome, region, sequence);
		fprintf(fa_fp, ">%s:%d-%d\n%s\n", region->chromosome, region->start, region->end, sequence->seq);

		fgets(buf, BUFFER, bed_fp);
	}

	fclose(fa_fp);
	fclose(bed_fp);


	return 0;
}


chromosomes_info_t * open_chromsize(const char * genome_dir) {
	FILE * fp;
	char path[BUFFER];
	strncpy(path, genome_dir, BUFFER);
	if(path[strlen(path)-1] != '/') {
		strcat(path, "/");
	}
	strcat(path, "chromsize");
	if(!(fp = fopen(path, "r"))) {
		fprintf(stdout, "[ERROR] can't open %s\n", path);
		return NULL;
	}
	char buf[BUFFER];
	char * token;

	fgets(buf, BUFFER, fp); //Read and ignore header
	//Count lines
	int dim = 0;
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		dim ++;
		fgets(buf, BUFFER, fp);
	}
	//printf("%d\n", dim);
	rewind(fp);

	char ** names;
	if( !(names=(char**)malloc(sizeof(char*)*dim)) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}
	int i;
	for(i=0; i<dim; i++) {
		if( !(names[i]=(char*)malloc(sizeof(char)*6)) ) {
			fprintf(stdout, "[ERROR] can't allocate\n");
			return NULL;
		}
	}

	int * sizes;
	if( !(sizes=malloc(sizeof(int)*dim)) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}

	fgets(buf, BUFFER, fp); //Read and ignore header
	i=0;
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, "\t");
		strncpy(names[i], token, 6);
		token = strtok(NULL, "\t");
		sizes[i] = atoi(token);
		i++;
		fgets(buf, BUFFER, fp);
	}

	fclose(fp);

	chromosomes_info_t * chrom_info = chromosome_info_init(names, sizes, dim);

	return chrom_info;
}
