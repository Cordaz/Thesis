#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "../utilities/argparse.h"
#include "../utilities/uthash.h"

typedef struct hashable_region_s {
	int start; // is the key
	int end;
	char strand;
	int score;
	UT_hash_handle hh;
} hashable_region_t;

void uthash_add_hashable_region(hashable_region_t ** t, int start, int end, char strand, int score) {
    hashable_region_t * h;

    h = malloc(sizeof(hashable_region_t));
    h->start = start;
	 h->end = end;
	 h->strand = strand;
	 h->score = score;
    HASH_ADD_INT( *t, start, h );
}

hashable_region_t * uthash_find_hashable_region(hashable_region_t ** t, int start) {
    hashable_region_t * h;

    HASH_FIND_INT( *t, &start, h );  /* s: output pointer */

	 return h;
}

void uthash_delete_all(hashable_region_t ** t) {
  hashable_region_t * cur, * tmp;

  HASH_ITER(hh, *t, cur, tmp) {
    HASH_DEL(*t, cur);  /* delete; users advances to next */
    free(cur);            /* optional- if you want to free  */
  }
}

void uthash_print_hashable_regions(hashable_region_t ** t, FILE * fp, char * chr) {
    hashable_region_t * h;

    for(h = *t; h != NULL; h=(hashable_region_t *)(h->hh.next)) {
        fprintf(fp, "%s\t%d\t%d\t%c\t%d\n", chr, h->start, h->end, h->strand, h->score);
    }
}

int uthash_start_sort(hashable_region_t * a, hashable_region_t * b) {
    return (a->start - b->start);
}

void uthash_sort_by_start(hashable_region_t ** t) {
    HASH_SORT(*t, uthash_start_sort);
}

#ifndef BUFFER
	#define BUFFER 512
#endif

static const char* const usage[] = {
	"faQuery [options] -f <in.fa> -i <in.kmers>",
	NULL
};

int main(int argc, const char * argv[]) {
	const char * fa_file = NULL;
	const char * kmers_file = NULL;
	int s = 8;

	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_STRING('i', "input", &kmers_file, "kmers input file"),
		OPT_STRING('f', "fasta", &fa_file, "fasta file to query"),
		OPT_INTEGER('s', "search-size", &s, "motif size to search (default 8)"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nSearches a fasta file for exact occurences of given list.", "\nFasta file must be in the form '>chr#:position' to retrieve the region\n");

	argc = argparse_parse(&argparse, argc, argv);

	if(fa_file == NULL || kmers_file == NULL) {
		fprintf(stdout, "[ERROR] must specify both fasta and input files\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	// General variables
	char buf[BUFFER+1];
	char buf2[BUFFER+1];
	int i;

	////////// Read input kmers
	// Open file
	FILE * kmers_fp;
	if( !(kmers_fp = fopen(kmers_file, "r")) ) {
		fprintf(stdout, "[ERROR] cannot open %s\n", kmers_file);
		return 1;
	}

	// Count and determine max length
	int max_kmers_len = 0;
	int len;
	int kmers_count = 0;
	fgets(buf, BUFFER, kmers_fp);
	while(!feof(kmers_fp)) {
		kmers_count++;
		len = strlen(buf) - 1; // Without '\n'
		if(len > max_kmers_len) {
			max_kmers_len = len;
		}
		fgets(buf, BUFFER, kmers_fp);
	}
	rewind(kmers_fp);

	// Alloc and read
	char ** kmers;
	if( !(kmers = (char**)malloc(sizeof(char*) * kmers_count)) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return 1;
	}
	int * kmers_len;
	if( !(kmers_len = (int *)malloc(sizeof(int) * kmers_count)) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return 1;
	}
	for(i=0; i<kmers_count; i++) {
		if( !(kmers[i] = (char*)malloc(sizeof(char) * max_kmers_len)) ) {
			fprintf(stdout, "[ERROR] cannot allocate\n");
			return 1;
		}
		fgets(buf, BUFFER, kmers_fp);
		kmers_len[i] = strlen(buf) - 1;
		strncpy(kmers[i], buf, kmers_len[i]);
		kmers[i][kmers_len[i]] = '\0';
		//printf("%s\t%d\n", kmers[i], kmers_len[i]);
	}
	fclose(kmers_fp);

	// Open fasta file
	FILE * fa_fp;
	if( !(fa_fp = fopen(fa_file, "r")) ) {
		fprintf(stdout, "[ERROR] cannot open %s\n", fa_file);
		return 1;
	}

	// Open bed output file
	FILE * bed_fp;
	char * extension;
	char bed_file[BUFFER + 1];
	strncpy(bed_file, kmers_file, BUFFER);
	extension = strrchr(bed_file, (int)'.');
	//printf("%s\n", extension);
	sprintf(extension, ".bed");
	//printf("%s\n", bed_file);
	if( !(bed_fp = fopen(bed_file, "w+")) ) {
		fprintf(stdout, "[ERROR] cannot open %s\n", bed_file);
		return 1;
	}

	// Ready hash_table
	hashable_region_t * region_table = NULL;
	hashable_region_t * region;
	/*
	uthash_add_hashable_region(&region_table, 1, 2, '+', 1);
	region = uthash_find_hashable_region(&region_table, 1);
	*/

	// Read fasta line by line
	char * ptr;
	char * token;
	char chr[6];
	char last_chr[6];
	int seq_start;
	int reg_start;
	fgets(buf, BUFFER, fa_fp); // Read header
	if(feof(fa_fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF in %s\n", fa_file);
		return 1;
	}
	fgets(buf2, BUFFER, fa_fp); // Read sequence
	strcpy(chr, "chr");
	for(i=3; buf[i+1] != ':'; i++) {
		chr[i] = buf[i+1];
	}
	strcpy(last_chr, chr);

	while(!feof(fa_fp)) {
		// Identify chromosome
		strcpy(chr, "chr");
		for(i=3; buf[i+1] != ':'; i++) {
			chr[i] = buf[i+1];
		}
		if(strcmp(last_chr, chr) != 0) {
			uthash_sort_by_start(&region_table);
			uthash_print_hashable_regions(&region_table, bed_fp, last_chr);
			uthash_delete_all(&region_table);
			strcpy(last_chr, chr);
		}
		token = NULL; // Used to just find seq_start only once per seq
		// For each query kmer
		for(i=0; i<kmers_count; i++) {
			ptr = strstr(buf2, kmers[i]);
			while(ptr) { // use if(ptr) and remove next ptr = strstr() call if one occurences per sequence maximum
				if(!token) {
					token = strtok(buf, ":");
					token = strtok(NULL, ":");
					//printf("%s\n", token);
					seq_start = atoi(token);
				}
				reg_start = seq_start + ( ptr - buf2 );
				//printf("%d\t%d\n", seq_start, reg_start);
				if( !(region = uthash_find_hashable_region(&region_table, reg_start)) ) {
					uthash_add_hashable_region(&region_table, reg_start, reg_start + s, '+', 1);
				} else {
					region->score += 1;
				}

				ptr = strstr(ptr+1, kmers[i]);
			}
		}

		fgets(buf, BUFFER, fa_fp); // Read header
		if(!feof(fa_fp)) {
			fgets(buf2, BUFFER, fa_fp); // Read sequence
		}
	}

	fclose(bed_fp);
	fclose(fa_fp);


	return 0;
}
