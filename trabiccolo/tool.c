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
#include "../utilities/FIFO.h"
#include "../utilities/argparse.h"

#ifdef SILENT
	#define printf(...)
#endif

#define BUFFER 256
#define MIN_ARGS 2
#define MAX_SUBS 2
#define BUILD 1
#define LOAD 2
#define FASTA 1
#define FASTQ 2

#define K 10

static const char* const usage[] = {
	"dsc [options] (-p path_pattern | -f format -i input -e experiment) -s search_size",
	NULL
};

graph_t * load_graph(const char *, int, int, int);
set_t * extend_right(set_t *, char *, int, int);

graph_t * build_graph(double, double, int);
int map_read(char *, int, int, graph_t *, fifo_t *);
int map_input_read(char *, int, int, graph_t *, fifo_t *);
node_t * get_successor(node_t *, int, char);

FILE * get_file_info(FILE *, int *, int *);

int main(int argc, const char * argv[]) {
	//Setup
	const char * pattern = NULL;
	const char * out_p = NULL;
	int input_format = 0;
	const char * input_file = NULL;
	const char * experiment_file = NULL;
	int s = 0;
	int k = K;
	int d = 0;
	int l;
	int num_of_subs = 2;
	int psm_arg = 0;
	int g = 0;
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
		OPT_INTEGER('k', "kmer", &k, "kmer length of graph (default 10)"),
		OPT_INTEGER('N', "num-subs", &num_of_subs, "Accepted substitution in approximate counting (default 2)"),
		OPT_BOOLEAN('d', "double", &d, "Consider double strand"),
		OPT_BOOLEAN('m', "psm", &psm_arg, "output Position Specific Matrix ext='.psm'"),
		OPT_STRING('n', "name", &out_p, "pattern to name output file (default 'pid.out')"),
		OPT_BOOLEAN('g', "graph", &g, "output graph (of experiment) ext='.graph'"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nSearch for approximate occurences of k-mer of length s in the adjDBG.\nEither specify the pattern of *.graph.* files (Required (1)) or specify Required (2) group", "\nNote that the grah is intended to be structured on k/2");

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
	int i, j, h;

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
	if(d) {
		fprintf(stdout, "                  @params d 'true'\n");
	} else {
		fprintf(stdout, "                  @params d 'false'\n");
	}
	fprintf(stdout, "                  @params k %d\n", k);
	fprintf(stdout, "                  @params N %d\n", num_of_subs);
	fprintf(stdout, "                  @output pattern '%s'\n", out_pattern);
	fprintf(stdout, "                  @output approximate count\n");
	if(psm_arg) {
		fprintf(stdout, "                  @output psm\n");
	}
	if(g) {
		fprintf(stdout, "                  @output graph\n");
	}


	int nodes = (int)pow((double)4, k/2);
	int edges = (int)pow((double)4, k/2+1);
	graph_t * dbg;

	//Either build or load the graph
	if(build_or_load == LOAD) {
		dbg = load_graph(pattern, k/2, nodes, edges);
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

		fifo_t * q;
		if ( !(q = init_queue(k/2)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}

		if( !(fp = fopen(experiment_file, "r")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", experiment_file);
			return 1;
		}

		int skip_line;
		if( !(fp = get_file_info(fp, &input_format, &l)) ) {
			return 1;
		}
		if (input_format == FASTA)
			skip_line = 2;
		else //is FASTQ
			skip_line = 4;

		char * read;
		if( !(read = (char*)malloc(sizeof(char) * (l+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}
		char * rev_read;
		if( !(rev_read = (char*)malloc(sizeof(char) * (l+1))) ) {
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
		fgets(buf, BUFFER, fp);
		i=0;
		while(!feof(fp)) {
			i++;
			if(i==2) {
				//This line has the read
				strncpy(read, buf, l); //Remove '\n', ensure length
				read[l] = '\0';
				//printf("%s\n", read);
				index = 0;
				while( (index = get_next_substring(read, index, k, &sublen)) != -1 ) {
					//Each substring should be mapped
					if (map_read(read+index, sublen, k, dbg, q)) {
						fprintf(stdout, "[ERROR] couldn't allocate\n");
						return 1;
					}
					index = index + sublen;
				}
			}
			if(i==skip_line) {
				i=0;
			}
			fgets(buf, BUFFER, fp);
		}

		fclose(fp);

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Processing of ChIP-seq complete\n");


		//// MAPPING CONTROL FILE
		if( !(fp = fopen(input_file, "r")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", input_file);
			return 1;
		}

		if( !(fp = get_file_info(fp, &input_format, &l)) ) {
			return 1;
		}
		if (input_format == FASTA)
			skip_line = 2;
		else //is FASTQ
			skip_line = 4;

		free(read);
		if( !(read = (char*)malloc(sizeof(char) * (l+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}
		free(rev_read);
		if( !(rev_read = (char*)malloc(sizeof(char) * (l+1))) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return 1;
		}

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Reading '%s'\n", input_file);

		//Read first line
		fgets(buf, BUFFER, fp);
		i=0;
		while(!feof(fp)) {
			i++;
			if(i==2) {
				//This line has the read
				strncpy(read, buf, l); //Remove '\n', ensure length
				read[l] = '\0';
				//printf("%s\n", read);
				index = 0;
				while( (index = get_next_substring(read, index, k, &sublen)) != -1 ) {
					//Each substring should be mapped
					if (map_input_read(read+index, sublen, k, dbg, q)) {
						fprintf(stdout, "[ERROR] couldn't allocate\n");
						return 1;
					}
					index = index + sublen;
				}
			}
			if(i==skip_line) {
				i=0;
			}
			fgets(buf, BUFFER, fp);
		}

		fclose(fp);

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Processing of Input complete\n");
	}

	// Getting total
	unsigned long total = 0;
	unsigned long total_input = 0;
	for(i=0; i<nodes; i++) {
		for(j=0; j<nodes; j++) {
			total += (unsigned long)dbg->nodes[i]->out_kstep[j]->count;
			total_input += (unsigned long)dbg->nodes[i]->out_kstep[j]->input_count;
		}
	}


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

	//Counters
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

	int dim = 1 + 3*s + 3*s/2 * 3*(s-1);
	set_t * subs;
	if( !(subs = initialize_set(dim, s)) ) {
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

	int flag;
	int flag2;

	//Approximate counting and psm
	for(i=0; i<expected_smer; i++) {
		clear(subs);
		//printf("%s\n", smers[i]);
		put(subs, smers[i]);
		if( !(subs = substitute(subs, smers[i], s, 0, num_of_subs)) ) {
			fprintf(stdout, "[ERROR] queue is not working\n");
			return 1;
		}
		flag = get(subs, smer);
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
					psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
				}

				flag2 = get(q, kmer);
			}

			//REVERSE
			if(d) {
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
						counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
						//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->count);
						input_counts[i] += (unsigned long)dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count;
						//printf("%d\n", dbg->nodes[hash_half0]->out_kstep[hash_half1]->input_count);

						for(j=0; j<s; j++) {
							psm[i][get_base_index(smer[j])][j] += dbg->nodes[hash_half0]->out_kstep[hash_half1]->count;
						}

						flag2 = get(q, kmer);
					}
				}
			}
			//END REVERSE

			flag = get(subs, smer);
		}

		//printf("%d\t%s\t%lu\t%lu\n", i, smers[i], counts[i], input_counts[i]);
	}


	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Occurences counted\n");

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

	//Output psm
	if(psm_arg) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: position specific matrix\n");
		strncpy(out_file, out_pattern, BUFFER);
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

	//Output graph
	if(g) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: graph\n");
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
		fprintf(fp_edges_out_k, "ID\tFrom\tTo\tCount\tInput_count\n");
		for(i=0; i<nodes; i++) {
			fprintf(fp_nodes, "%d\t%s\n", i, dbg->nodes[i]->seq);

			for(j=0; j<4; j++) {
				fprintf(fp_edges_out, "%d\t%s\t%s\n", dbg->nodes[i]->out[j]->id, dbg->nodes[i]->seq, dbg->nodes[i]->out[j]->to->seq);
			}

			for(j=0; j<nodes; j++) {
				fprintf(fp_edges_out_k, "%d\t%s\t%s\t%d\t%d\n", dbg->nodes[i]->out_kstep[j]->id, dbg->nodes[i]->seq, dbg->nodes[i]->out_kstep[j]->to->seq, dbg->nodes[i]->out_kstep[j]->count, dbg->nodes[i]->out_kstep[j]->input_count);
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

int map_read(char * read, int l, int k, graph_t * dbg, fifo_t * q) {
	int i;
	node_t * n = dbg->nodes[hash(read, k/2)]; //Get starting node
	node_t * n0;
	q = enqueue(q, n);
	//printf("%x enqueued\n", n);
	for (i=1; i<(k/2); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1)); //Get the node corresponding to the right overlapping kmer
		q = enqueue(q, n);
		//printf("%x enqueued\n", n);
	}
	//Now start to add edges for contiguos kmer
	for(; i<(l-k/2+1); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1));
		n0 = dequeue(q);
		//printf("n0: %s, n: %s, read: %s\n", n0->seq, n->seq, read);
		dbg->nodes[n0->id]->out_kstep[n->id]->count += 1;

		q = enqueue(q, n);
	}

	return 0;
}

int map_input_read(char * read, int l, int k, graph_t * dbg, fifo_t * q) {
	int i;
	node_t * n = dbg->nodes[hash(read, k/2)]; //Get starting node
	node_t * n0;
	q = enqueue(q, n);
	//printf("%x enqueued\n", n);
	for (i=1; i<(k/2); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1)); //Get the node corresponding to the right overlapping kmer
		q = enqueue(q, n);
		//printf("%x enqueued\n", n);
	}
	//Now start to add edges for contiguos kmer
	for(; i<(l-k/2+1); i++) {
		n = get_successor(n, k/2, *(read+i+k/2-1));
		n0 = dequeue(q);
		//printf("n0: %s, n: %s, read: %s\n", n0->seq, n->seq, read);
		dbg->nodes[n0->id]->out_kstep[n->id]->input_count += 1;


		q = enqueue(q, n);
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

FILE * get_file_info(FILE * fp, int * format, int * reads_len) {
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

	fgets(buf, BUFFER, fp);
	//printf("%s\n", buf);
	*reads_len = strlen(buf) - 1; //Ignoring '\n'

	rewind(fp);

	return fp;
}
