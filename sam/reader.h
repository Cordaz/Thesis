#ifndef SAM_READER_H
	#define SAM_READER_H

	#include <stdio.h>
	#include <string.h>
	#include <malloc.h>
	#include <stdlib.h>

	typedef struct sam_header_s {
		char chromosomes[24][3];
		int chromosomes_size[24];
	} sam_header_t;

	typedef struct sam_file_s {
		FILE * fp;
		sam_header_t * header;
	} sam_file_t;

	typedef struct sam_line_s {
		char chromosome[3];
		int flag;
		int pos;
		int read_len;
	} sam_line_t;

	sam_file_t * sam_open(char *, char *, int);

	sam_line_t * sam_readline(sam_file_t *, sam_line_t *, char *, int, int);

#endif
