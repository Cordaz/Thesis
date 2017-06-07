#ifndef GENOME_H
	#define GENOME_H

	/*
	 * struct of info of the genome
	 * @field dim number of chromosome
	 * @field names chromosome names (chr*)
	 * @field sizes length of chromosomes
	 *
	 */
	typedef struct chromosomes_info_s {
		int dim;
		char ** names;
		int * sizes;
	} chromosomes_info_t;

	/*
	 * struct of chromosome
	 * @field name (chr*)
	 * @field index (as in chromosome_info struct or in the BAM/SAM header)
	 * @field seq sequence
	 * @field size length
	 *
	 */
	typedef struct chromosome_s {
		char name[6];
		int index;
		char * seq;
		int size;
	} chromosome_t;

	/*
	 * struct genome
	 * @field path path/to/genome/dir/
	 * @field chromosome (active)
	 * @field info chromosome info struct
	 *
	 */
	typedef struct genome_s {
		char path[512+1];
		chromosome_t * chromosome;
		chromosomes_info_t * info;
	} genome_t;

	/*
	 * struct for a region
	 * @field chromosome name
	 * @field chrom_index as defined in chrom_info
	 * @field start starting position
	 * @field end end position
	 * @field strand '+','-','.'
	 *
	 */
	typedef struct region_s {
		char chromosome[6];
		int chrom_index;
		int start;
		int end;
		char strand;
	} region_t;

	/*
	 * Initialize the genome struct (see above)
	 * @param genome_path path/to/genome/dir/ (containing chr*.fa files)
	 * @param chrom_info struct (already initialized)
	 * @return genome
	 *
	 */
	genome_t * genome_init(char *, chromosomes_info_t *);

	/*
	 * Clear the genome struct
	 * @param geome struct to be cleared
	 *
	 */
	void genome_clear(genome_t *);

	/*
	 * load a chromosome in memory, storing in genome struct its pointer
	 * @param genome struct
	 * @param chromosome_name
	 * @param chrom_index (as defined above)
	 * @return genome with a new active chromosome
	 *
	 */
	genome_t * genome_load_chromosome(genome_t *, char *, int);

	/*
	 * get the index corresponding to the provided chromosome name according to chrom_info struct
	 * @param chrom_info struct
	 * @param chromosome_name
	 * @return index, -1 if not found
	 *
	 */
	int get_chrom_index(chromosomes_info_t *, char *);

	/*
	 * Initialize chrom_info struct from a BAM/SAM header (decomposed in its part) or from separated info (i.e. read by a file)
	 * @param names array of chromosome names
	 * @param sizes array of chromosome sizes
	 * @param dim dimension of the genome (i.e. num of chromosome)
	 * @return chrom_info struct initialized, define chromosome index
	 *
	 */
	chromosomes_info_t * chromosome_info_init(char **, int *, int);

	/*
	 * Vlear chromosome info struct
	 * @param chrom_info struct to be cleared
	 *
	 */
	void chromosome_info_clear(chromosomes_info_t *);

#endif
