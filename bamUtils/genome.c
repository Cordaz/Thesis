#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <ctype.h>

#include "genome.h"

genome_t * genome_init(char * path, chromosomes_info_t * info) {
	if(path[strlen(path)-1] != '/') path[strlen(path)-1] = '/';

	genome_t * genome;

	if( !(genome = (genome_t *)malloc(sizeof(genome_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}

	strcpy(genome->path, path);
	genome->chromosome = NULL;

	genome->info = info;

	return genome;
}

void genome_clear(genome_t * genome) {
	free(genome->chromosome->seq);
	free(genome->chromosome);
	free(genome);
}

genome_t * genome_load_chromosome(genome_t * genome, char * chromosome_name) {
	if(!genome->chromosome) {
		if(!(genome->chromosome = (chromosome_t*)malloc(sizeof(chromosome_t)))) {
			fprintf(stdout, "[ERROR] can't allocate\n");
			return NULL;
		}
		genome->chromosome->seq = NULL;
	}

	if(genome->chromosome->seq) {
		free(genome->chromosome->seq);
	}

	int index;
	index = get_chrom_index(genome->info, chromosome_name);
	genome->chromosome->index = index;
	int size;
	size = genome->info->sizes[index];
	genome->chromosome->size = size;
	strcpy(genome->chromosome->name, chromosome_name);

	if( !(genome->chromosome->seq = (char*)malloc(sizeof(char) * (size + 1))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}

	char buf[512];
	strcpy(buf, genome->path);
	strcat(buf, chromosome_name);
	strcat(buf, ".fa");

	FILE * fp;
	if( !(fp=fopen(buf, "r")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", buf);
		return NULL;
	}

	char c;

	int i=0;
	fgets(buf, 512, fp); //Read and ignore first line
	while(!feof(fp)) {
		fread(&c, 1, 1, fp);
		if(c != '\n') {
			if(islower(c)) c = toupper(c);
			genome->chromosome->seq[i] = c;
			i++;
		}
	}


	fclose(fp);

	return genome;

}

int get_chrom_index(chromosomes_info_t * chrom_info, char * chrom_name) {
	int i;
	for(i=0; i<chrom_info->dim; i++) {
		if( strcmp(chrom_info->names[i], chrom_name) == 0 ) return i;
	}
	return -1;
}


chromosomes_info_t * chromosome_info_init(char ** names, int * sizes, int dim) {
	int i;
	chromosomes_info_t * chrom_info;

	if( !(chrom_info = (chromosomes_info_t*)malloc(sizeof(chromosomes_info_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}
	if( !(chrom_info->names = (char **)malloc(sizeof(char*) * dim)) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}
	for(i=0; i<dim; i++) {
		if( !(chrom_info->names[i] = (char *)malloc(sizeof(char) * 6)) ) {
			fprintf(stdout, "[ERROR] can't allocate\n");
			return NULL;
		}
	}

	if( !(chrom_info->sizes = (int *)malloc(sizeof(int) * dim)) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}

	for(i=0; i<dim; i++) {
		strncpy(chrom_info->names[i], names[i], 5);
		chrom_info->sizes[i] = sizes[i];
	}
	chrom_info->dim = dim;

	return chrom_info;
}


void chromosome_info_clear(chromosomes_info_t * chrom_info) {
	int i;
	for(i=0; i<chrom_info->dim; i++) {
		free(chrom_info->names[i]);
	}
	free(chrom_info->names);
	free(chrom_info->sizes);
	free(chrom_info);
}
