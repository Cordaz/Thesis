#ifndef GENOME_H
	#define GENOME_H

	typedef struct chromosomes_info_s {
		int dim;
		char ** names;
		int * sizes;
	} chromosomes_info_t;

	typedef struct chromosome_s {
		char name[6];
		int index;
		char * seq;
		int size;
	} chromosome_t;

	typedef struct genome_s {
		char path[512+1];
		chromosome_t * chromosome;
		chromosomes_info_t * info;
	} genome_t;

	typedef struct region_s {
		char chromosome[6];
		int start;
		int end;
	} region_t;

	/*
	 * Path to genome/ expected
	 *
	 */
	genome_t * genome_init(char *, chromosomes_info_t *);

	void genome_clear(genome_t *);

	genome_t * genome_load_chromosome(genome_t *, char *);

	int get_chrom_index(chromosomes_info_t *, char *);

	chromosomes_info_t * chromosome_info_init(char **, int *, int);

	void chromosome_info_clear(chromosomes_info_t *);

#endif
