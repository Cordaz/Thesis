#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "../utilities/data_structures.h"
#include "../utilities/my_lib.h"
#include "../utilities/set.h"
#include "../utilities/argparse.h"

#ifdef SILENT
	#define printf(...)
#endif

#define BUFFER 256
#define MIN_ARGS 2
#define MAX_SUBS 2

static const char* const usage[] = {
	"dsc [options] -p path_pattern -s search_size",
	NULL
};

graph_t * load_graph(const char *, int, int, int);
set_t * extend_right(set_t *, char *, int, int);

int main(int argc, const char * argv[]) {
	//Setup
	const char * pattern = NULL;
	const char * out_p = NULL;
	int s = 0;
	int k = 10;
	int psm_arg = 0;
	int c_arg = 0;
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Mandatory"),
		OPT_STRING('p', "pattern", &pattern, "pattern to .graph.*, .psm and .approx.cnt files"),
		OPT_INTEGER('s', "search", &s, "size of motifs to retrieve (k/2 < s < k)"),
		OPT_GROUP("Optional"),
		OPT_INTEGER('k', "kmer", &k, "kmer length of graph (default 10)"),
		OPT_BOOLEAN('c', "count", &c_arg, "output the approximate count (default not done)"),
		OPT_BOOLEAN('m', "psm", &psm_arg, "output Position Specific Matrix ext='.psm'"),
		OPT_STRING('n', "name", &out_p, "pattern to name output file (default 'pid.out')"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nSearch for approximate occurences of k-mer of length s in the adjDBG.", "\nNote that the grah is intended to be structured on k/2");
	if(argc < MIN_ARGS) {
		argparse_usage(&argparse);
		return 0;
	}
	argc = argparse_parse(&argparse, argc, argv);

	if(!s) {
		fprintf(stdout, "[ERROR] search size must be specified\n\n");
		argparse_usage(&argparse);
		return 0;
	}
	if( s <= k/2 || s > k ) {
		fprintf(stdout, "[ERROR] search size must be in k/2 < s <= k\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(!pattern) {
		fprintf(stdout, "[ERROR] pattern must be specified\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	int pid = getpid();

	char out_pattern[BUFFER+1];
	if(!out_p) {
		sprintf(out_pattern, "%d.out", pid);
	} else {
		strcpy(out_pattern, out_p);
	}
	time_t rawtime;
	struct tm * timeinfo;
	//End setup

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Starting with PID %d\n", pid);
	fprintf(stdout, "                  @params s %d\n", s);
	fprintf(stdout, "                  @params pattern '%s'\n", pattern);
	fprintf(stdout, "                  @params k %d\n", k);
	if(c_arg || psm_arg) { //OR of each output
		fprintf(stdout, "                  @output pattern '%s'\n", out_pattern);
	}
	if(c_arg) {
		fprintf(stdout, "                  @output approximate count\n");
	}
	if(psm_arg) {
		fprintf(stdout, "                  @output psm\n");
	}


	int nodes = (int)pow((double)4, k/2);
	int edges = (int)pow((double)4, k/2+1);
	graph_t * dbg = load_graph(pattern, k/2, nodes, edges);
	if(!dbg) {
		return 1;
	}
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Loaded graph in memory\n");

	int i, j, h;

	//estimate number of s-mer
	int expected_smer;
	expected_smer = (int)pow((double)4, s);

	//// Matrix for position specific count
	int *** psm; //Position Specific Matrix
	if( !(psm = (int***)malloc(sizeof(int**) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_smer; i++) {
		if( !(psm[i] = (int**)malloc(sizeof(int*) * 4)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		for(j=0; j<4; j++) {
			if( !(psm[i][j] = (int*)malloc(sizeof(int) * s)) ) {
				fprintf(stdout, "[ERROR] couldn't allocate\n");
				return 1;
			}
			for(h=0; h<s; h++) {
				psm[i][j][h] = 0;
			}
		}
	}

	unsigned long * counts;
	if( !(counts = (unsigned long*)malloc(sizeof(unsigned long) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	unsigned long * input_counts;
	if( !(input_counts = (unsigned long*)malloc(sizeof(unsigned long) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_smer; i++) {
		counts[i] = 1;
		input_counts[i] = 1;
		//Pseudo-count
	}

	int q_len = (int)pow( (double)4, k-s );
	set_t * q;
	if( !(q = initialize_set(q_len, k)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	int dim = 1 + 3*s;
	set_t * one;
	if( !(one = initialize_set(dim, s)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	dim = dim + 3*s/2 * 3*(s-1);
	set_t * two;
	if( !(two = initialize_set(dim, s)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	unsigned long total = expected_smer; //Fills the gap of Pseudo-count
	unsigned long total_input = expected_smer;

	char ** smers;
	if( !(smers = (char**)malloc(sizeof(char*) * (expected_smer) )) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_smer; i++) {
		if( !(smers[i] = (char*)malloc(sizeof(char) * (s+1) )) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		rev_hash(i, s, smers[i]);
		//printf("%s\n", smers[i]);
	}

	char * kmer;
	if( !(kmer = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	char * smer;
	if( !(smer = (char*)malloc(sizeof(char)*(s+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	char * rev_smer;
	if( !(rev_smer = (char*)malloc(sizeof(char)*(s+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	char * half0, * half1;
	if( !(half0 = (char*)malloc(sizeof(char) * (k/2+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	if( !(half1 = (char*)malloc(sizeof(char) * (k/2+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	int hash_half0, hash_half1;

	int p;
	int flag;
	int flag2;

	for(i=0; i<expected_smer; i++) {
		clear(one);
		clear(two);
		//printf("%s\n", smers[i]);
		if( !(one = substitute_one(one, smers[i], s)) ) {
			fprintf(stdout, "[ERROR] queue is not working\n");
			return 1;
		}
		flag = get(one, smer);
		while( flag != -1 ) {
			//printf("\t%s\n", smer);
			two = substitute_one(two, smer, s);
			/*
			str_dequeue(two, kmer);
			while( strcmp(kmer, "") ) {
				printf("\t\t%s\n", kmer);
				str_dequeue(two, kmer);
			}
			*/
			flag = get(one, smer);
		}

		flag = get(two, smer);
		while( flag != -1 ) {
			//printf("\t%s\n", smer);
			clear(q);

			if( !(q = extend_right(q, smer, k-s, k)) ) {
				fprintf(stdout, "[ERROR] couldn't allocate\n");
				return 1;
			}
			flag2 = get(q, kmer);
			while( flag2 != -1 ) {
				//printf("%s\t", kmer);
				strncpy(half0, kmer, k/2);
				strncpy(half1, kmer+k/2, k/2);
				hash_half0 = hash(half0, k/2);
				hash_half1 = hash(half1, k/2);
				//printf("%s\t%s\t", half0, half1);
				counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
				//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->count);
				input_counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
				//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count);

				for(j=0; j<s; j++) {
					psm[i][get_base_index(kmer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
				}

				flag2 = get(q, kmer);
			}

			//Reverse complementary TODO map directly on the graph
			reverse_kmer(smer, rev_smer, s);
			//printf("%s\t%s\n", smer, rev_smer);

			if(!is_palyndrome(smer, rev_smer)) {
				clear(q);

				if( !(q = extend_right(q, rev_smer, k-s, k)) ) {
					fprintf(stdout, "[ERROR] couldn't allocate\n");
					return 1;
				}
				flag2 = get(q, kmer);
				while( flag2 != -1 ) {
					//printf("%s\t", kmer);
					strncpy(half0, kmer, k/2);
					strncpy(half1, kmer+k/2, k/2);
					hash_half0 = hash(half0, k/2);
					hash_half1 = hash(half1, k/2);
					//printf("%s\t%s\t", half0, half1);
					counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
					//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->count);
					input_counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
					//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count);

					for(j=0; j<s; j++) {
						psm[i][get_base_index(kmer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
					}

					flag2 = get(q, kmer);
				}
			}


			flag = get(two, smer);
		}

		//printf("\n");

		//printf("%d\t%s\t%lu\t%lu\n", i, smers[i], counts[i], input_counts[i]);
	}


	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Occurences counted\n");

	//// OUTPUT

	char buf[BUFFER+1];
	char out_file[BUFFER+1];
	FILE * fp;
	if(c_arg) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: count\n");
		strncpy(out_file, out_pattern, BUFFER);
		sprintf(buf, ".s%d", s);
		strcat(out_file, buf);
		strcat(out_file, ".approx.cnt");
		if( !(fp = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		double freq;
		double freq_input;
		double diff;
		double diff_log2;

		fprintf(fp, "k-mer\tIP_count\tIP_freq\tInput_count\tInput_freq\tdiff\tdiff_log2\n");
		for(i=0; i<expected_smer; i++) {
			rev_hash(i, s, kmer);
			freq = (double)counts[i]/(double)total;
			freq_input = (double)input_counts[i]/(double)total_input;
			diff = freq/freq_input;
			diff_log2 = log2(diff);
			fprintf(fp, "%s\t%lu\t%lf\t%lu\t%lf\t%lf\t%lf\n", kmer, counts[i], freq, input_counts[i], freq_input, diff, diff_log2);
		}
		fclose(fp);
	}

	if(psm_arg) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: position specific matrix\n");
		strncpy(out_file, out_pattern, BUFFER);
		sprintf(buf, ".s%d", s);
		strcat(out_file, buf);
		strcat(out_file, ".psm");
		if( !(fp = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}

		for(i=0; i<expected_smer; i++) {
			fprintf(fp, "%s\n", smers[i]);
			fprintf(fp, "A");
			for(j=0; j<s; j++) {
				fprintf(fp, "\t%d", psm[i][0][j]);
			}
			fprintf(fp, "\nC");
			for(j=0; j<s; j++) {
				fprintf(fp, "\t%d", psm[i][1][j]);
			}
			fprintf(fp, "\nG");
			for(j=0; j<s; j++) {
				fprintf(fp, "\t%d", psm[i][2][j]);
			}
			fprintf(fp, "\nT");
			for(j=0; j<s; j++) {
				fprintf(fp, "\t%d", psm[i][3][j]);
			}
			fprintf(fp, "\n\n");
		}

		fclose(fp);
	}



	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Completed\n");


	return 0;
}


graph_t * load_graph(const char * pattern, int k, int nodes, int edges) {
	int i;
	char buf[BUFFER+1];
	char * token;
	const char sep[2] = "\t";
	char args[5][BUFFER+1];
	//Allocate graph
	graph_t * dbg;
	if( !(dbg = (graph_t*)malloc(sizeof(graph_t))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}
	if( !(dbg->nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	if( !(dbg->edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	FILE * fp;
	char file_path[BUFFER+1];
	node_t * n;
	int node_id;
	strcpy(file_path, pattern);
	strcat(file_path, ".graph.nodes");
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Load nodes
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<2 ) { //Expecting two args
			strcpy(args[i], token);

			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 2 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0'; //Remove '\n'
		node_id = atoi(args[0]);
		if( !(n = create_node(node_id, args[1], nodes)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		dbg->nodes[node_id] = n;

		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	edge_t * e;
	int id;
	int n0_hash;
	int n1_hash;

	//Load one-step out edges
	strcpy(file_path, pattern);
	strcat(file_path, ".graph.edges.out");
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<3 ) { //Expecting three args
			strcpy(args[i], token);
			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 3 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0';
		args[2][k] = '\0';
		id = atoi(args[0]);
		n0_hash = hash(args[1], k);
		n1_hash = hash(args[2], k);
		if( !( e = create_edge( dbg->nodes[n0_hash], dbg->nodes[n1_hash], id ) ) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		dbg->edges[id] = e;
		dbg->nodes[n0_hash]->out[get_base_index(args[2][k-1])] = e;
		dbg->nodes[n1_hash]->in[get_base_index(args[1][0])] = e;

		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	//Load k/2-step out edges
	strcpy(file_path, pattern);
	sprintf(buf, ".graph.edges.out.%dstep", k);
	strcat(file_path, buf);
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<5 ) { //Expecting five args
			strcpy(args[i], token);
			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 5 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0';
		args[2][k] = '\0';
		id = atoi(args[0]);
		n0_hash = hash(args[1], k);
		n1_hash = hash(args[2], k);
		if( !( e = create_edge( dbg->nodes[n0_hash], dbg->nodes[n1_hash], id ) ) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		e->count = atoi(args[3]);
		e->input_count = atoi(args[4]);
		dbg->nodes[n0_hash]->out_kstep[n1_hash] = e;
		dbg->nodes[n1_hash]->in_kstep[n0_hash] = e;
		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	return dbg;
}

set_t * extend_right(set_t * q, char * str, int times, int l) {
	if (times == 0) {
		put(q, str);
		return q;
	}
	if(times < 0)
		return q;

	char * support;
	if( !(support = (char*)malloc(sizeof(char) * (l+1))) )
		return NULL;
	strcpy(support, str);
	int i = 0;
	while(str[i] != '\0')
		i++;
	support[i+1] = '\0';
	support[i] = 'A';
	q = extend_right(q, support, times-1, l);
	support[i] = 'C';
	q = extend_right(q, support, times-1, l);
	support[i] = 'G';
	q = extend_right(q, support, times-1, l);
	support[i] = 'T';
	q = extend_right(q, support, times-1, l);

	free(support);

	return q;
}
