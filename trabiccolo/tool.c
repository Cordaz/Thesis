/*
@author Andrea Corneo
@year 2017

Main part of 'trabiccolo' project. Still not released.

*/

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "../utilities/data_structures.h"
#include "../utilities/my_lib.h"
#include "../utilities/string_FIFO.h"
#include "../utilities/node_FIFO.h"
#include "../utilities/argparse.h"
#include "../utilities/list.h"
#include "../utilities/sorting.h"
#include "../utilities/bit_array.h"
#include "../utilities/distance.h"

#ifdef SILENT
	#define printf(...)
#endif

#ifndef BUFFER
	#define BUFFER 512
#endif
#ifndef SEQUENCE_BUFFER
	#define SEQUENCE_BUFFER 1048576
#endif
#define MIN_ARGS 2
#define MAX_SUBS 2
#define BUILD 1
#define LOAD 2
#define FASTA 1
#define FASTQ 2

#define K 10

static const char* const usage[] = {
	"tool [options] (-p path_pattern | -i input -e experiment) -s search_size",
	NULL
};

const char BASES[4] = {'A','C','G','T'};

int ** update_psm(int **, int, char *, int);

graph_t * load_graph(const char *, int, int, int, unsigned long *, unsigned long *, unsigned long *, unsigned long *);
string_FIFO_t * extend_right(string_FIFO_t *, char *, int, int);

graph_t * build_graph(double, double, int);
int map_read(char *, int, int, graph_t *, node_FIFO_t *, int *, int, unsigned long *);
int map_input_read(char *, int, int, graph_t *, node_FIFO_t *, int *, int, unsigned long *);
node_t * get_successor(node_t *, int, char);

FILE * get_file_info(FILE *, int *);

int main(int argc, const char * argv[]) {
	char rscript_path[BUFFER];
	strcpy(rscript_path, argv[0]);

	//Setup
	const char * pattern = NULL;
	const char * out_p = NULL;
	int input_format = 0;
	const char * input_file = NULL;
	const char * experiment_file = NULL;
	const char * adapters_file = NULL;
	int s = 0;
	int k = K;
	int S = 0;
	int r_arg = 0;
	int l;
	int num_of_subs = 1;
	int psm_arg = 0;
	int g_arg = 0;
	int build_or_load = 0;
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Mandatory"),
		OPT_INTEGER('s', "search", &s, "size of motifs to retrieve (k/2 < s < k)"),
		OPT_GROUP("Required (1)"),
		OPT_STRING('p', "pattern", &pattern, "pattern to .graph.*, .psm and .approx.cnt files"),
		OPT_GROUP("Required (2)"),
		OPT_STRING('i', "input", &input_file, "Input file (Control)"),
		OPT_STRING('e', "experiment", &experiment_file, "Experiment file"),
		OPT_GROUP("Optional"),
		OPT_BOOLEAN('R', "region-count", &r_arg, "Count occurences in region instead of total occurences"),
		OPT_STRING('a', "adapters", &adapters_file, "Adapters file"),
		OPT_INTEGER('k', "kmer", &k, "kmer length of graph (default 10)"),
		OPT_INTEGER('N', "num-subs", &num_of_subs, "Accepted substitution in approximate counting (default 2, maximum 2)"),
		OPT_BOOLEAN('S', "single", &S, "Consider single strand only"),
		OPT_BOOLEAN('m', "psm", &psm_arg, "output Position Specific Matrix ext='.psm'"),
		OPT_STRING('n', "name", &out_p, "pattern to name output file (default 'pid.out')"),
		OPT_BOOLEAN('g', "graph", &g_arg, "output graph (of experiment) ext='.graph'"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nSearch for approximate occurences of k-mer of length s in the adjDBG.\nEither specify the pattern of *.graph.* files (Required (1)) or specify Required (2) group", "\nNote that the graph is intended to be structured on k/2");

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

	if(num_of_subs > 2) {
		fprintf(stdout, "[ERROR] maximum number of accepted substitution must be <= 2\n\n");
		argparse_usage(&argparse);
		return 0;
	}

	//printf("%s\t%s\t%s\n", format, input_file, experiment_file);
	if(!pattern) {
		if(!experiment_file) {
			fprintf(stdout, "[ERROR] required group 1 or 2 should be specified\n\n");
			argparse_usage(&argparse);
			return 0;
		}
	}
	if(!experiment_file) {
		if(!pattern) {
			fprintf(stdout, "[ERROR] required group 1 or 2 should be specified\n\n");
			argparse_usage(&argparse);
			return 0;
		}
	}

	if(pattern) {
		build_or_load = LOAD;
	} else if (input_file && experiment_file) {
		build_or_load = BUILD;
	} else {
		fprintf(stdout, "[ERROR] required group 1 or 2 should be specified\n\n");
		argparse_usage(&argparse);
		return 0;
	}


	int pid = getpid();

	FILE * fp;
	char buf[BUFFER+1];
	int i, j, h, g;

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
	if(build_or_load == LOAD)
		fprintf(stdout, "                  @params pattern '%s'\n", pattern);
	if(build_or_load == BUILD) {
		fprintf(stdout, "                  @params input '%s'\n", input_file);
		fprintf(stdout, "                  @params experiment '%s'\n", experiment_file);
	}
	if(S) {
		fprintf(stdout, "                  @params S 'true'\n");
	} else {
		fprintf(stdout, "                  @params S 'false'\n");
	}
	if(r_arg) {
		fprintf(stdout, "                  @params R 'true': counting occurences in region\n");
	} else {
		fprintf(stdout, "                  @params R 'false': counting total occurences of kmer\n");
	}
	fprintf(stdout, "                  @params k %d\n", k);
	fprintf(stdout, "                  @params N %d\n", num_of_subs);
	fprintf(stdout, "                  @output pattern '%s'\n", out_pattern);
	fprintf(stdout, "                  @output approximate count\n");
	if(psm_arg) {
		fprintf(stdout, "                  @output psm\n");
	}
	if(g_arg) {
		fprintf(stdout, "                  @output graph\n");
	}


	int nodes = (int)pow((double)4, k/2);
	int edges = (int)pow((double)4, k/2+1);
	graph_t * dbg;

	unsigned long kmer_total = 0;
	unsigned long kmer_total_input = 0;
	unsigned long reads_total = 0;
	unsigned long reads_total_input = 0;


	//Either build or load the graph
	if(build_or_load == LOAD) {
		dbg = load_graph(pattern, k/2, nodes, edges, &reads_total, &reads_total_input, &kmer_total, &kmer_total_input);
		if(!dbg) {
			return 1;
		}
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Loaded graph in memory\n");
	} else if (build_or_load == BUILD) {
		dbg = build_graph(nodes, edges, k/2);
		if (!dbg) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Empty De Bruijn graph built\n");
		fprintf(stdout, "                  created %d nodes\n", (int)nodes);
		fprintf(stdout, "                  created %d 1-step edges\n", (int)edges);
		fprintf(stdout, "                  created %d %d-step edges\n", (int)(nodes*nodes), k/2);

		int * last_seen;
		if( !(last_seen = (int*)malloc(sizeof(int)*nodes*nodes)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		for(i=0; i<nodes*nodes; i++) {
			last_seen[i] = -1;
		}

		int read_index;

		node_FIFO_t * q;
		if ( !(q = init_queue(k/2)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}

		if( !(fp = fopen(experiment_file, "r")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", experiment_file);
			return 1;
		}

		int skip_line;
		if( !(fp = get_file_info(fp, &input_format)) ) {
			return 1;
		}
		if (input_format == FASTA)
			skip_line = 2;
		else //is FASTQ
			skip_line = 4;

		char * read;
		if( !(read = (char*)malloc(sizeof(char) * (SEQUENCE_BUFFER+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}
		char * seq_buf;
		if( !(seq_buf = (char*)malloc(sizeof(char) * (SEQUENCE_BUFFER+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}

		int index;
		int sublen;

		//printf("%d\t%d\n", l, input_format);

		//// OPENING READS FILE
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Reading '%s'\n", experiment_file);


		//Read first line
		fgets(seq_buf, SEQUENCE_BUFFER, fp);
		i=0;
		read_index = 0;
		while(!feof(fp)) {
			i++;
			if(i==2) {
				//This line has the read
				strncpy(read, seq_buf, SEQUENCE_BUFFER+1);
				l = strlen(read) - 1;  //Remove '\n', ensure length
				read[l] = '\0';
				//printf("%s\n", read);
				index = 0;
				while( (index = get_next_substring(read, index, k, &sublen)) != -1 ) {
					//Each substring should be mapped
					if (map_read(read+index, sublen, k, dbg, q, last_seen, read_index, &kmer_total)) {
						fprintf(stdout, "[ERROR] couldn't allocate\n");
						return 1;
					}
					index = index + sublen;
				}
				read_index++;
			}
			if(i==skip_line) {
				i=0;
			}
			fgets(seq_buf, SEQUENCE_BUFFER, fp);
		}

		fclose(fp);

		reads_total = read_index;

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Processing of ChIP-seq complete\n");


		//// MAPPING CONTROL FILE
		if( !(fp = fopen(input_file, "r")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", input_file);
			return 1;
		}

		if( !(fp = get_file_info(fp, &input_format)) ) {
			return 1;
		}
		if (input_format == FASTA)
			skip_line = 2;
		else //is FASTQ
			skip_line = 4;

		free(read);
		if( !(read = (char*)malloc(sizeof(char) * (SEQUENCE_BUFFER+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Reading '%s'\n", input_file);

		//Read first line
		fgets(seq_buf, SEQUENCE_BUFFER, fp);
		read_index = 0;
		for(i=0; i<nodes*nodes; i++) {
			last_seen[i] = -1;
		}
		i=0;

		while(!feof(fp)) {
			i++;
			if(i==2) {
				//This line has the read
				strncpy(read, seq_buf, SEQUENCE_BUFFER+1);
				l = strlen(read) - 1;  //Remove '\n', ensure length
				read[l] = '\0';
				//printf("%s\n", read);
				index = 0;
				while( (index = get_next_substring(read, index, k, &sublen)) != -1 ) {
					//Each substring should be mapped
					if (map_input_read(read+index, sublen, k, dbg, q, last_seen, read_index, &kmer_total_input)) {
						fprintf(stdout, "[ERROR] couldn't allocate\n");
						return 1;
					}
					index = index + sublen;
				}
				read_index++;
			}
			if(i==skip_line) {
				i=0;
			}
			fgets(seq_buf, SEQUENCE_BUFFER, fp);
		}

		fclose(fp);

		reads_total_input = read_index;

		free(last_seen);
		free(read);
		free(seq_buf);

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Processing of Input complete\n");
	}

	char * kmer;
	if( !(kmer = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	//Reading adapters
	list_t * adapters = NULL;

	if(adapters_file) {
		if( !(fp = fopen(adapters_file, "r")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", adapters_file);
			return 1;
		}

		fgets(buf, BUFFER, fp);
		i = 0;
		while(!feof(fp)) {
			i++;
			if(i==2) {
				i=0;
				l = strlen(buf) - 1;
				//printf("%s\n", buf);

				for(j=0; j<l-s+1; j++) {
					strncpy(kmer, buf+j, s);
					kmer[s] = '\0';
					if(!search(adapters, kmer)) {
						adapters = add(adapters, kmer);
						if(!(adapters) ) {
							fprintf(stdout, "[ERROR] couldn't allocate memory\n");
							return 1;
						}
					}
				}
			}
			fgets(buf, BUFFER, fp);
		}

		/*
		list_t * t;
		t = adapters;
		while(t) {
			printf("%s\n", t->str);
			t = t->next;
		}
		*/

		fclose(fp);
	}

	unsigned long total;
	unsigned long total_input;

	const int R_arg = r_arg;
	if(R_arg) {
		total = reads_total;
		total_input = reads_total_input;
	} else {
		total = kmer_total;
		total_input = kmer_total_input;
	}

	unsigned long count;
	unsigned long input_count;
	double freq;
	double freq_input;
	double diff;
	double diff_log2;


	//estimate number of s-mer
	int expected_smer;
	expected_smer = (int)pow((double)4, s);

	/*
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
	*/

	//Counters
	unsigned long ** counts;
	if( !(counts = (unsigned long**)malloc(sizeof(unsigned long*) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	unsigned long ** input_counts;
	if( !(input_counts = (unsigned long**)malloc(sizeof(unsigned long*) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_smer; i++) {
		if( !(counts[i] = (unsigned long*)malloc(sizeof(unsigned long) * (num_of_subs + 1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		if( !(input_counts[i] = (unsigned long*)malloc(sizeof(unsigned long) * (num_of_subs + 1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		for(g=0; g<=num_of_subs; g++) {
			counts[i][g] = 0;
			input_counts[i][g] = 1;
		}
	}

	int q_len = (int)pow( (double)4, k-s );
	string_FIFO_t * q;
	if( !(q = initialize_set(q_len, k)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	//int dim = 1 + 3*s + 3*s/2 * 3*(s-1);
	int dim = 3*s * 3*s;
	string_FIFO_t * subs;
	if( !(subs = initialize_set(dim, s)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	string_FIFO_t * subs_reserve;
	if( !(subs_reserve = initialize_set(dim, s)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

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

	int flag;
	int flag2;

	//Approximate counting and psm
	for(i=0; i<expected_smer; i++) {
		clear(subs_reserve);
		put(subs_reserve, smers[i]);
		for(g=0; g <= num_of_subs; g++) {
			clear(subs);
			flag = get(subs_reserve, smer);
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
					if(R_arg) {
						count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
						input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
					} else {
						count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
						input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_input_count;
					}

					freq = (double)count/total;
					freq_input = (double)(input_count)/total_input;
					if(freq >= freq_input) {
						counts[i][g] += count;
						input_counts[i][g] += input_count;
						/*
						// This part is used to create psm matrix for each entry.
						// Substituted by a "greedy" matrix composition of top motifs
						for(j=0; j<s; j++) {
							if(R_arg) {
								psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
							} else {
								psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
							}
						}
						*/
					}
					flag2 = get(q, kmer);
				}

				//REVERSE
				if(!S) {
					reverse_kmer(smer, rev_smer, s);
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
							if(R_arg) {
								count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
								input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
							} else {
								count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
								input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_input_count;
							}

							freq = (double)count/total;
							freq_input = (double)(input_count)/total_input;
							if(freq >= freq_input) {
								counts[i][g] += count;
								input_counts[i][g] += input_count;
								/*
								for(j=0; j<s; j++) {
									if(R_arg) {
										psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
									} else {
										psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
									}
								}
								*/
							}

							flag2 = get(q, kmer);
						}
					}
				}
				//END REVERSE

				if(g != num_of_subs) {
					if( !(subs = substitute(subs, smer, s, 0, 1)) ) {
						fprintf(stdout, "[ERROR] queue is not working\n");
						return 1;
					}
				}

				flag = get(subs_reserve, smer);
			}

			clear(subs_reserve);
			flag = get(subs, smer);
			while( flag != -1 ) {
				put(subs_reserve, smer);
				flag = get(subs, smer);
			}
		}

	}


	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Occurences counted\n");

	double * measures;
	if( !(measures = (double*)malloc(sizeof(double) * expected_smer)) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return 1;
	}

	//// OUTPUT
	char out_file[BUFFER+1];

	//Output approximate count
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Generating output: count\n");
	strncpy(out_file, out_pattern, BUFFER);
	strcat(out_file, ".approx.cnt");
	if( !(fp = fopen(out_file, "w+")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", out_file);
		return 1;
	}

	fprintf(fp, "k-mer\tIP_count_0\tIP_freq_0\tInput_count_0\tInput_freq_0\tdiff_0\tdiff_log2_0");
	if(num_of_subs >= 1) {
		fprintf(fp, "\tIP_count_1\tIP_freq_1\tInput_count_1\tInput_freq_1\tdiff_1\tdiff_log2_1");
	}
	if(num_of_subs >= 2) {
		fprintf(fp, "\tIP_count_2\tIP_freq_2\tInput_count_2\tInput_freq_2\tdiff_2\tdiff_log2_2");
	}
	if(num_of_subs >= 1) {
		fprintf(fp, "\tsum_of_freq\tsum_of_freq_input\tlog2_ratio_freq");
	}
	fprintf(fp, "\n");
	for(i=0; i<expected_smer; i++) {
		smer = smers[i];
		/*
		if( search(adapters, smer) ) {
			for(j=0; j<s; j++) {
				smer[j] = tolower(smer[j]);
			}
		}
		*/
		fprintf(fp, "%s", smer);
		// Printing no sub
		freq = (double)counts[i][0]/(double)total;
		freq_input = (double)input_counts[i][0]/(double)total_input;
		diff = freq / freq_input;
		diff_log2 = log2( diff );
		measures[i] = diff_log2;

		fprintf(fp, "\t%lu\t%lf\t%lu\t%lf\t%lf\t%lf", counts[i][g], freq, input_counts[i][g], freq_input, diff, diff_log2);
		// Printing subs
		for(g=1; g <= num_of_subs; g++) {
			freq = (double)counts[i][g]/(double)total;
			freq_input = (double)input_counts[i][g]/(double)total_input;
			diff = freq / freq_input;
			diff_log2 = log2( diff );

			fprintf(fp, "\t%lu\t%lf\t%lu\t%lf\t%lf\t%lf", counts[i][g], freq, input_counts[i][g], freq_input, diff, diff_log2);
		}

		if(num_of_subs >= 1) {
			fprintf(fp, "\t%lf", measures[i]);
		}

		fprintf(fp, "\n");

	}
	fclose(fp);


	/***************** GREEDY PSM GENERATION *****************************************************************************/

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Generating greedy matrix\n");

	int * sorted = sort_index_decreasing(measures, expected_smer);

	//// Matrix for position specific count
	int to_select = 10; // TODO as arg
	int *** psm; //Position Specific Matrix
	if( !(psm = (int***)malloc(sizeof(int**) * to_select)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<to_select; i++) {
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

	int * motifs;
	if( !(motifs = (int*)malloc(sizeof(int) * to_select)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	int * positions;
	if( !(positions = (int*)malloc(sizeof(int) * to_select)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	int limit = num_of_subs + 1;

	bit_array_t * flagged = init_bit_array(expected_smer);
	// Check for adapters only if existing list
	if(adapters) {
		for(i=0; i<expected_smer; i++) {
			if( search(adapters, smers[sorted[i]]) ) {
				if(!(flagged = set_on_bit(flagged, i))) {return 1;}
				for(j=0; j<expected_smer; j++) {
					if(dist(smers[sorted[i]], smers[sorted[j]], s, limit) < limit) {
						if(!(flagged = set_on_bit(flagged, j))) {return 1;}
					}
				}
			}
		}
	}

	int processed = 0;
	int selected = -1;
	int hash_val;
	int rev_index;

	int * minus, * plus;
	if( !(minus = (int*)malloc(sizeof(int) * to_select)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	if( !(plus = (int*)malloc(sizeof(int) * to_select)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	while (processed < to_select) {
		do {
			selected++;
			//printf("\t%s\t%d\t%d\n", smers[sorted[selected]], selected, sorted[selected]);
		} while(selected < expected_smer && get_bit(flagged, selected));
		if (selected >= expected_smer) {
			to_select = processed; // If less are computed then ajust output to this
			break;
		}
		motifs[processed] = sorted[selected];
		positions[processed] = selected;
		//printf("\t%d\n", motifs[processed]);
		psm[processed] = update_psm(psm[processed], s, smers[sorted[selected]], counts[sorted[selected]][0]);
		plus[processed] = 1;
		minus[processed] = 0;
		flagged = set_on_bit(flagged, selected);
		if(!flagged) { return 1; }

		// Flag also the reverse comp
		reverse_kmer(smers[sorted[selected]], rev_smer, s);
		hash_val = hash(rev_smer, s);
		rev_index = get_pos(sorted, selected+1, expected_smer, hash_val);
		//printf("%s\t%s\t%d\t%d\n", smers[sorted[selected]], smers[hash_val], hash_val, rev_index);
		if(rev_index >= 0 && rev_index < expected_smer) {
			flagged = set_on_bit(flagged, rev_index);
			if(!flagged) { return 1; }
		}

		for(i=selected+1; i<expected_smer; i++) {
			if (dist(smers[sorted[selected]], smers[sorted[i]], s, limit) < limit) {
				if (measures[sorted[i]] > 0) { // FIXME measures o freq/freq_input con 0 sostituzioni?
					psm[processed] = update_psm(psm[processed], s, smers[sorted[i]], counts[sorted[i]][0]);
					plus[processed]++;
				} else {
					minus[processed]++;
				}
				flagged = set_on_bit(flagged, i);
				if(!flagged) { return 1; }
				// Also here flag rev comp
				reverse_kmer(smers[sorted[i]], rev_smer, s);
				hash_val = hash(rev_smer, s);
				rev_index = get_pos(sorted, i+1, expected_smer, hash_val);
				//printf("%s\t%s\t%d\t%d\n", smers[sorted[selected]], smers[hash_val], hash_val, rev_index);
				if(rev_index >= 0 && rev_index < expected_smer) {
					flagged = set_on_bit(flagged, rev_index);
					if(!flagged) { return 1; }
				}
			}
		}

		processed++;
	}
	/////////////////////// END NEW
	//Output psm
	if(psm_arg) {
		//printf("%s\n", rscript_path);
		char * str_ptr = strrchr(rscript_path, '/');
		strcpy(str_ptr, "/");
		//printf("%s\n", str_ptr);

		char rscript[BUFFER];
		sprintf(rscript, "Rscript %slogo.R ", rscript_path);
		char cmd[BUFFER];

		// Creating subdirectory
		struct stat st = {0};

		if (stat(out_pattern, &st) == -1) {
			 mkdir(out_pattern, 0700);
		}

		char dir[BUFFER];
		snprintf(dir, BUFFER, "%s/", out_pattern);

		char html_file[BUFFER];
		snprintf(html_file, BUFFER, "%s.motifs.html", out_pattern);

		char kmer_file[BUFFER];
		FILE * kmer_fp;

		FILE * html_fp;
		if( !(html_fp = fopen(html_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", html_file);
			return 1;
		}
		fprintf(html_fp, "<html>\n<head>\n<title>%s - Motifs</title>\n</head>\n\n<body>\n", out_pattern);

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: position specific matrix\n");

		for(i=0; i<to_select; i++) {
			strncpy(kmer_file, dir, BUFFER);
			sprintf(buf, "%d_%s.kmers", i+1, smers[motifs[i]]);
			strcat(kmer_file, buf);
			if( !(kmer_fp = fopen(kmer_file, "w+")) ) {
				fprintf(stdout, "[ERROR] can't open %s\n", kmer_file);
				return 1;
			}


			// Print list of used kmer
			clear(subs_reserve);
			put(subs_reserve, smers[motifs[i]]);
			for(g=0; g <= num_of_subs; g++) {
				clear(subs);
				flag = get(subs_reserve, smer);
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
						if(R_arg) {
							count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
							input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
						} else {
							count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
							input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_input_count;
						}

						freq = (double)count/total;
						freq_input = (double)(input_count)/total_input;
						if(freq >= freq_input) {
							fprintf(kmer_fp, "%s\n", kmer);
						}
						flag2 = get(q, kmer);
					}

					//REVERSE
					if(!S) {
						reverse_kmer(smer, rev_smer, s);
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
								if(R_arg) {
									count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
									input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
								} else {
									count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_count;
									input_count = (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->kmer_input_count;
								}

								freq = (double)count/total;
								freq_input = (double)(input_count)/total_input;
								if(freq >= freq_input) {
									fprintf(kmer_fp, "%s\n", kmer);
								}

								flag2 = get(q, kmer);
							}
						}
					}
					//END REVERSE

					if(g != num_of_subs) {
						if( !(subs = substitute(subs, smer, s, 0, 1)) ) {
							fprintf(stdout, "[ERROR] queue is not working\n");
							return 1;
						}
					}

					flag = get(subs_reserve, smer);
				}

				clear(subs_reserve);
				flag = get(subs, smer);
				while( flag != -1 ) {
					put(subs_reserve, smer);
					flag = get(subs, smer);
				}
			}




			strncpy(out_file, dir, BUFFER);
			sprintf(buf, "%d_%s.psm", i+1, smers[motifs[i]]);
			strcat(out_file, buf);
			if( !(fp = fopen(out_file, "w+")) ) {
				fprintf(stdout, "[ERROR] can't open %s\n", out_file);
				return 1;
			}
			fprintf(html_fp, "<h3>%s</h3>\n", smers[motifs[i]]);
			fprintf(html_fp, "<p>Starting kmer index: %d</p>\n", positions[i]+1);
			fprintf(html_fp, "<p>Used: %d kmers<br/>Discarded: %d kmers</p>\n", plus[i], minus[i]);
			fprintf(html_fp, "<table>");
			for(g=1; g<=num_of_subs; g++) {
				fprintf(html_fp, "<tr><td>IP_count_%d</td><td>IP_freq_%d</td><td>Input_count_%d</td><td>Input_freq_%d</td></tr>", g, g, g, g);
				freq = (double)counts[motifs[i]][g]/(double)total;
				freq_input = (double)input_counts[motifs[i]][g]/(double)total_input;
				fprintf(html_fp, "<tr><td>%lu</td><td>%lf</td><td>%lu</td><td>%lf</td></tr>", counts[motifs[i]][g], freq, input_counts[motifs[i]][g], freq_input);
			}
			fprintf(html_fp, "<tr><td colspan='2'>log2_ratio_freq</td><td colspan='2'>%lf</td></tr>", measures[motifs[i]]);
			fprintf(html_fp, "</table>");
			fprintf(html_fp, "<table>\n");

			for(h=0; h<4; h++) {
				fprintf(html_fp, "<tr><td><b>%c</b></td>", BASES[h]);
				for(j=0; j<s-1; j++) {
					fprintf(fp, "%d\t", psm[i][h][j]);
					fprintf(html_fp, "<td>%d</td>", psm[i][h][j]);
				}
				fprintf(fp, "%d\n", psm[i][h][s-1]);
				fprintf(html_fp, "<td>%d</td></tr>\n", psm[i][h][s-1]);
			}
			fclose(fp);

			fprintf(html_fp, "</table>\n");

			// Run R script to plot sequence logo
			strcpy(cmd, rscript);
			strcat(cmd, dir);
			strcat(cmd, buf);
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
			fprintf(stdout, "Exectuing '%s'\n", cmd);
			system(cmd);

			fprintf(html_fp, "<img src='%s%s.png' width='800px' height='300px'/>\n\n", dir, buf);
		}

		fprintf(html_fp, "</body>\n</html>");
		fclose(html_fp);
	}

	//Output graph
	if(g_arg) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: graph\n");
		strncpy(out_file, out_pattern, BUFFER);
		strcat(out_file, ".graph.info");
		FILE * fp_info;
		if( !(fp_info = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		fprintf(fp_info, "@reads_total\t%lu\n", reads_total);
		fprintf(fp_info, "@reads_total_input\t%lu\n", reads_total_input);
		fprintf(fp_info, "@kmer_total\t%lu\n", kmer_total);
		fprintf(fp_info, "@kmer_total_input\t%lu\n", kmer_total_input);
		fclose(fp_info);
		strncpy(out_file, out_pattern, BUFFER);
		strcat(out_file, ".graph.nodes");
		FILE * fp_nodes;
		if( !(fp_nodes = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		strncpy(out_file, out_pattern, BUFFER);
		strcat(out_file, ".graph.edges.out");
		FILE * fp_edges_out;
		if( !(fp_edges_out = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		strncpy(out_file, out_pattern, BUFFER);
		strcat(out_file, ".graph.edges.out");
		sprintf(buf, ".%dstep", k/2);
		strcat(out_file, buf);
		FILE * fp_edges_out_k;
		if( !(fp_edges_out_k = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		fprintf(fp_nodes, "ID\tseq\n");
		fprintf(fp_edges_out, "ID\tFrom\tTo\n");
		fprintf(fp_edges_out_k, "ID\tFrom\tTo\tCount\tKmer_count\tInput_count\tKmer_input_count\n");
		for(i=0; i<nodes; i++) {
			fprintf(fp_nodes, "%d\t%s\n", i, dbg->nodes[i]->seq);

			for(j=0; j<4; j++) {
				fprintf(fp_edges_out, "%d\t%s\t%s\n", dbg->nodes[i]->out[j]->id, dbg->nodes[i]->seq, dbg->nodes[i]->out[j]->to->seq);
			}

			for(j=0; j<nodes; j++) {
				fprintf(fp_edges_out_k, "%d\t%s\t%s\t%d\t%d\t%d\t%d\n", dbg->nodes[i]->out_kstep[j]->id, dbg->nodes[i]->seq, dbg->nodes[i]->out_kstep[j]->to->seq, dbg->nodes[i]->out_kstep[j]->count, dbg->nodes[i]->out_kstep[j]->kmer_count, dbg->nodes[i]->out_kstep[j]->input_count, dbg->nodes[i]->out_kstep[j]->kmer_input_count);
			}
		}
		fclose(fp_nodes);
		fclose(fp_edges_out);
		fclose(fp_edges_out_k);
	}

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Completed\n");

	free(smer);


	return 0;
}



int ** update_psm(int ** psm, int s, char * smer, int count) {
	int i;
	for (i=0; i<s; i++) {
		psm[get_base_index(smer[i])][i] += count;
	}

	return psm;
}



graph_t * load_graph(const char * pattern, int k, int nodes, int edges, unsigned long * reads_total, unsigned long * reads_total_input, unsigned long * kmer_total, unsigned long * kmer_total_input) {
	int i;
	char buf[BUFFER+1];
	char * token;
	const char sep[2] = "\t";
	char args[7][BUFFER+1];

	FILE * fp;
	char file_path[BUFFER+1];
	//Getting infos
	strncpy(file_path, pattern, BUFFER);
	strcat(file_path, ".graph.info");
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	fgets(buf, BUFFER, fp);
	if(feof(fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF in %s\n", file_path);
	}
	char * endptr;
	token = strtok(buf, sep); //read first arg - ignore
	token = strtok(NULL, sep); //second arg
	*reads_total = strtoul(token, &endptr, 10);
	fgets(buf, BUFFER, fp);
	if(feof(fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF in %s\n", file_path);
	}
	token = strtok(buf, sep); //read first arg - ignore
	token = strtok(NULL, sep); //second arg
	*reads_total_input = strtoul(token, &endptr, 10);
	fgets(buf, BUFFER, fp);
	if(feof(fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF in %s\n", file_path);
	}
	token = strtok(buf, sep); //read first arg - ignore
	token = strtok(NULL, sep); //second arg
	*kmer_total = strtoul(token, &endptr, 10);
	fgets(buf, BUFFER, fp);
	if(feof(fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF in %s\n", file_path);
	}
	token = strtok(buf, sep); //read first arg - ignore
	token = strtok(NULL, sep); //second arg
	*kmer_total_input = strtoul(token, &endptr, 10);


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


	node_t * n;
	int node_id;
	strncpy(file_path, pattern, BUFFER);
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
	strncpy(file_path, pattern, BUFFER);
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
	strncpy(file_path, pattern, BUFFER);
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
		while( token && i<7 ) { //Expecting seven args
			strcpy(args[i], token);
			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 7 arguments, more found. Exit.\n");
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
		e->kmer_count = atoi(args[4]);
		e->input_count = atoi(args[5]);
		e->kmer_input_count = atoi(args[6]);
		dbg->nodes[n0_hash]->out_kstep[n1_hash] = e;
		dbg->nodes[n1_hash]->in_kstep[n0_hash] = e;
		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	return dbg;
}

string_FIFO_t * extend_right(string_FIFO_t * q, char * str, int times, int l) {
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

int map_read(char * read, int l, int k, graph_t * dbg, node_FIFO_t * q, int * last_seen, int read_index, unsigned long * kmer_total) {
	int i;
	int fullhash;
	node_t * n = dbg->nodes[hash(read, k/2)]; //Get starting node
	node_t * n0;
	q = enqueue(q, n);
	//printf("%x enqueued\n", n);
	for (i=1; i<(k/2); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1)); //Get the node corresponding to the right overlapping kmer
		q = enqueue(q, n);
		*kmer_total = *kmer_total + 1;
		//printf("%x enqueued\n", n);
	}
	//Now start to add edges for contiguos kmer
	for(; i<(l-k/2+1); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1));
		n0 = dequeue(q);
		//printf("n0: %s, n: %s, read: %s\n", n0->seq, n->seq, read);
		fullhash = (n0->id << k) + n->id;
		if(read_index > last_seen[fullhash]) {
			dbg->nodes[n0->id]->out_kstep[n->id]->count += 1;
			last_seen[fullhash] = read_index;
		}
		dbg->nodes[n0->id]->out_kstep[n->id]->kmer_count += 1;

		q = enqueue(q, n);
		*kmer_total = *kmer_total + 1;
	}

	return 0;
}

int map_input_read(char * read, int l, int k, graph_t * dbg, node_FIFO_t * q, int * last_seen, int read_index, unsigned long * kmer_total_input) {
	int i;
	int fullhash;
	node_t * n = dbg->nodes[hash(read, k/2)]; //Get starting node
	node_t * n0;
	q = enqueue(q, n);
	//printf("%x enqueued\n", n);
	for (i=1; i<(k/2); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1)); //Get the node corresponding to the right overlapping kmer
		q = enqueue(q, n);
		//printf("%x enqueued\n", n);
		*kmer_total_input = *kmer_total_input + 1;
	}
	//Now start to add edges for contiguos kmer
	for(; i<(l-k/2+1); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1));
		n0 = dequeue(q);
		//printf("n0: %s, n: %s, read: %s\n", n0->seq, n->seq, read);
		fullhash = (n0->id << k) + n->id;
		if(read_index > last_seen[fullhash]) {
			dbg->nodes[n0->id]->out_kstep[n->id]->input_count += 1;
			last_seen[fullhash] = read_index;
		}
		dbg->nodes[n0->id]->out_kstep[n->id]->kmer_input_count += 1;

		q = enqueue(q, n);
		*kmer_total_input = *kmer_total_input + 1;
	}

	return 0;
}

node_t * get_successor(node_t * n, int k, char c) {
	int i;
	for(i=0; i<4; i++) {
		if( n->out[i]->to->seq[k-1] == c )
			return n->out[i]->to;
	}
	return NULL;
}

graph_t * build_graph(double nodes, double edges, int k) {
	int i, j;
	int mask_hash = (int)nodes - 1;

	graph_t * dbg;
	if ( !(dbg = (graph_t*)malloc(sizeof(graph_t))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}
	if ( !(dbg->nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}
	if ( !(dbg->edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	char ** kmer;
	if( !(kmer = (char**)malloc(sizeof(char*)*nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}
	for(i=0; i<nodes; i++) {
		if ( !(kmer[i] = (char*)malloc(sizeof(char) * (k+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
	}
	char * support;
	if ( !(support = (char*)malloc(sizeof(char) * (2*k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	node_t * n;

	//Create nodes
	for (i=0; i<nodes; i++) {
		rev_hash(i, k, kmer[i]);
		if( !(n = create_node(i, kmer[i], nodes)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		dbg->nodes[i] = n;
	}

	int n0_hash, n1_hash;
	edge_t * e;

	//Create 1-step edges
	for (i=0; i<edges; i++) {
		n0_hash = i >> 2;
		n1_hash = i & mask_hash;

		if( !(e = create_edge(dbg->nodes[n0_hash], dbg->nodes[n1_hash], i)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}

		dbg->edges[i] = e;
		dbg->nodes[n0_hash]->out[(i%4)] = e;
		dbg->nodes[n1_hash]->in[n0_hash >> 2*(k-1)] = e;

	}

	//Create k/2-step edges
	for (i=0; i<nodes; i++) {
		for(j=0; j<nodes; j++) {
			strcpy(support, dbg->nodes[i]->seq);
			strcat(support, dbg->nodes[j]->seq);
			if( !(e = create_edge(dbg->nodes[i], dbg->nodes[j], hash(support, k*2) )) ) {
				fprintf(stdout, "[ERROR] couldn't allocate memory\n");
				return NULL;
			}
			dbg->nodes[i]->out_kstep[j] = e;
			dbg->nodes[j]->in_kstep[i] = e;
		}
	}

	for(i=0; i<nodes; i++) {
		free(kmer[i]);
	}
	free(kmer);

	return dbg;
}

FILE * get_file_info(FILE * fp, int * format) {
	char buf[BUFFER+1];
	fgets(buf, BUFFER, fp);
	//printf("%s\n", buf);
	//First line should start with '@' for fastq, '>' for fasta
	if(buf[0] == '@')
		*format = FASTQ;
	else if(buf[0] == '>')
		*format = FASTA;
	else {
		fprintf(stdout, "[ERROR] format not recognized\n");
		return NULL;
	}

	//Now get reads length
	if(feof(fp)) {
		fprintf(stdout, "[ERROR] unexpected EOF\n");
		return NULL;
	}

	rewind(fp);

	return fp;
}
