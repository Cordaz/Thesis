#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>

#include "reader.h"

/*
int main (int argc, char * argv[]) {
	char buffer[512];
	sam_file_t * samfile = sam_open(argv[1], buffer, 512);
	if(!samfile) return 1;
	sam_line_t * samline = NULL;
	samline = sam_readline(samfile, samline, buffer, 512, 100);
	printf("%s\t%d\n", samfile->header->chromosomes[1], samfile->header->chromosomes_size[1]);
	printf("%d\t%s\t%d\t%d\n", samline->flag, samline->chromosome, samline->pos, samline->read_len);
}
*/


sam_file_t * sam_open(char * path, char * buf, int buf_len) {
	FILE * fp;
	if( !(fp = fopen(path, "r")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", path);
		return NULL;
	}

	sam_file_t * samfile;
	if( !(samfile = (sam_file_t*)malloc(sizeof(sam_file_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}
	if( !(samfile->header = (sam_header_t*)malloc(sizeof(sam_header_t))) ) {
		fprintf(stdout, "[ERROR] can't allocate\n");
		return NULL;
	}
	samfile->fp = fp;

	char * token;

	fpos_t position;

	const char tab_sep[2] = "\t";
	const char colon_sep[2] = ":";

	fgetpos(fp, &position);
	fgets(buf, buf_len, fp);
	int i = 0;
	while (1) {
		if(buf[0] != '@') break;

		token = strtok(buf, tab_sep);
		if(strcmp(token, "@SQ") == 0) {
			token = strtok(NULL, tab_sep);
			if(token[0] == 'S' && token[1] == 'N') {
				strncpy(samfile->header->chromosomes[i], token+6, 3); //Skip SN:chr

				token = strtok(NULL, tab_sep);
				if(token[0] == 'L' && token[1] == 'N') {
					samfile->header->chromosomes_size[i] = atoi(token+3);
					i++;
				}
			}
		}
		fgetpos(fp, &position);
		fgets(buf, buf_len, fp);
	}

	//First line read, seek back
	fsetpos(fp, &position);

	return samfile;
}

sam_line_t * sam_readline(sam_file_t * samfile, sam_line_t * samline, char * buf, int buf_len, int extension) {
	if(!samline) {
		if( !(samline = (sam_line_t*)malloc(sizeof(sam_line_t))) ) {
			fprintf(stdout, "[ERROR] can't allocate\n");
			return NULL;
		}
	}

	char * token;
	const char tab_sep[2] = "\t";
	const char colon_sep[2] = ":";

	fgets(buf, buf_len, samfile->fp);
	if(feof(samfile->fp)) return NULL;

	token = strtok(buf, tab_sep);
	//Ignore first field
	token = strtok(NULL, tab_sep);
	samline->flag = atoi(token);
	token = strtok(NULL, tab_sep);
	strncpy(samline->chromosome, token, 3);
	token = strtok(NULL, tab_sep);
	samline->pos = atoi(token);
	if(!extension) {
		//Ignore 5 field
		token = strtok(NULL, tab_sep);
		token = strtok(NULL, tab_sep);
		token = strtok(NULL, tab_sep);
		token = strtok(NULL, tab_sep);
		token = strtok(NULL, tab_sep);

		token = strtok(NULL, tab_sep);
		samline->read_len = strlen(token);
	} else {
		samline->read_len = extension;
	}


	return samline;
}
