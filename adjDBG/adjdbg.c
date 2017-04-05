#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "../utilities/data_structures.h"
#include "../utilities/my_lib.h"
#include "../utilities/FIFO.h"
#include "../utilities/argparse.h"

#ifdef SILENT
	#define printf(...)
#endif

#define IN_FORMAT 1
#define K_ARG 2
#define K 10
#define L 34
#define L_ARG 4
#define O_ARG 6
#define C_FILE 7
#define IN_FILE 8
#define MIN_ARGS 3

#define BUFFER 256

#define FASTA 1
#define FASTQ 2

#define MAX_SUBS 2

// ARGUMENTS
static const char* const usage[] = {
	"adjdbg [options] -f fatsa|fastq -i input_file -e experiment_file",
	NULL
};


///////// FUNCTIONS PROTOTYPES

graph_t * build_graph(double, double, int);
int map_read(char *, int, int, graph_t *, fifo_t *);
int map_input_read(char *, int, int, graph_t *, fifo_t *);
node_t * get_successor(node_t *, int, char);


int main (int argc, const char * argv[]) {

	int k = K;
	int l = L;
	char * format = NULL;
	int input_format = 0;
	const char * input_file = NULL;
	const char * control_file = NULL;
	int c = 0;
	char out_file[BUFFER+1];
	int psm_arg = 0;
	int g = 0;

	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Mandatory"),
		OPT_STRING('f', "format", &format, "format of files to read"),
		OPT_STRING('i', "input", &control_file, "Input file (Control)"),
		OPT_STRING('e', "experiment", &input_file, "Experiment file"),
		OPT_GROUP("Optional"),
		OPT_INTEGER('k', "kmer", &k, "kmer length (default 10)"),
		OPT_INTEGER('l', "length", &l, "reads length (default 34)"),
		OPT_BOOLEAN('c', "count", &c, "output approximate count ext='.approx.cnt'"),
		OPT_BOOLEAN('m', "psm", &psm_arg, "output Position Specific Matrix ext='.psm'"),
		OPT_BOOLEAN('g', "graph", &g, "output graph (of experiment) ext='.graph'"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nBuilds a De Bruijn graph on k/2-mer, then maps input and experiment on the graph. Calculate approximate count and Position Specific Matrix.", "\nOutput files are names as experiment_file.k#.specific_extension");
	if(argc < MIN_ARGS) {
		argparse_usage(&argparse);
		return 0;
	}
	argc = argparse_parse(&argparse, argc, argv);

	if(strcmp(format, "fasta") == 0)
		input_format = FASTA;
	else if(strcmp(format, "fastq") == 0)
		input_format = FASTQ;

	if(!input_format) {
		fprintf(stdout, "[ERROR] format is mandatory\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(!input_file || !control_file) {
		fprintf(stdout, "[ERROR] -i input_file and -e experiment_file should be specified\n");
		argparse_usage(&argparse);
		return 0;
	}

	int pid = getpid();
	time_t rawtime;
	struct tm * timeinfo;



	//// -----------------------
	int i, j;
	int index;
	int sublen;
	char buf[BUFFER+1];
	char * read;
	if( !(read = (char*)malloc(sizeof(char) * (l+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return 1;
	}
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Starting with PID %d\n", pid);
	fprintf(stdout, "                  @params k %d\n", k);
	fprintf(stdout, "                  @params l %d\n", l);
	fprintf(stdout, "                  @params format %s\n", format);
	fprintf(stdout, "                  @params input %s\n", control_file);
	fprintf(stdout, "                  @params experiment %s\n", input_file);
	if(c) {
		fprintf(stdout, "                  @output approximate count\n");
	}
	if(psm_arg) {
		fprintf(stdout, "                  @output psm\n");
	}
	if(g) {
		fprintf(stdout, "                  @output graph\n");
	}

	//// BUILD EMPTY GRAPH
	double nodes = pow((double)4, (double)(k/2));
	double edges = pow((double)4, (double)((k/2)+1));
	graph_t * dbg = build_graph(nodes, edges, k/2);
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

	//// INIT QUEUE
	fifo_t * q;
	if ( !(q = init_queue(k/2)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}

	FILE * fp;


	//// OPENING READS FILE
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Reading %s\n", input_file);

	if( !(fp = fopen(input_file, "r")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", input_file);
		return 1;
	}
	int skip_line;
	if (input_format == FASTA)
		skip_line = 2;
	else //is FASTQ
		skip_line = 4;


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
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Reading %s\n", control_file);

	if( !(fp = fopen(control_file, "r")) ) {
		fprintf(stdout, "[ERROR] can't open %s\n", control_file);
		return 1;
	}
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


	//// WORKING ON FILLED GRAPH
	int expected_sub = 3*k/2;
	char ** substituted;
	if( !(substituted = (char**)malloc(sizeof(char*)*expected_sub)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_sub; i++) {
		if( !(substituted[i] = (char*)malloc(sizeof(char) * (k/2 + 1) )) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
	}
	char ** substituted_adj;
	if( !(substituted_adj = (char**)malloc(sizeof(char*)*expected_sub)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<expected_sub; i++) {
		if( !(substituted_adj[i] = (char*)malloc(sizeof(char) * (k/2 + 1) )) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
	}

	char * kmer;
	char * revkmer;
	if( !(revkmer = (char*)malloc(sizeof(char)*(k/2+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	char * adjkmer;
	char * adjrevkmer;
	if( !(adjrevkmer = (char*)malloc(sizeof(char)*(k/2+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	char * reference;
	if( !(reference = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	char * revreference;
	if( !(revreference = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	char * support;
	if( !(support = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	char * revsupport;
	if( !(revsupport = (char*)malloc(sizeof(char)*(k+1))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	node_t * substituted_node;
	node_t * rev_substituted_node;
	int h, m, x;
	unsigned long * counts;
	if( !(counts = (unsigned long*)malloc(sizeof(unsigned long) * nodes * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	unsigned long * counts_input;
	if( !(counts_input = (unsigned long*)malloc(sizeof(unsigned long) * nodes * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	unsigned long count;
	unsigned long count_input;
	unsigned long total = 0;
	unsigned long total_input = 0;

	int subhash;
	int revsubhash;
	int adjhash;
	int revadjhash;
	int fullhash;
	int revhash;


	//// Matrix for position specific count
	int *** psm; //Position Specific Matrix
	if( !(psm = (int***)malloc(sizeof(int**) * nodes * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate\n");
		return 1;
	}
	for(i=0; i<nodes*nodes; i++) {
		if( !(psm[i] = (int**)malloc(sizeof(int*) * 4)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate\n");
			return 1;
		}
		for(j=0; j<4; j++) {
			if( !(psm[i][j] = (int*)malloc(sizeof(int) * k)) ) {
				fprintf(stdout, "[ERROR] couldn't allocate\n");
				return 1;
			}
			for(h=0; h<k; h++) {
				psm[i][j][h] = 0;
			}
		}
	}

	//// Approximate count of k-mer
	for(i=0; i<nodes; i++) {
		kmer = dbg->nodes[i]->seq;
		substitute_all(kmer, substituted, k/2);
		reverse_kmer(kmer, revkmer, k/2);

		for(m=0; m<nodes; m++) {
			strcpy(reference, kmer);
			adjkmer = dbg->nodes[i]->out_kstep[m]->to->seq;
			strcat(reference, adjkmer); //Reference k-mer
			count = (unsigned long)1; //Pseudo-count
			count_input = (unsigned long)1; //Pseudo-count
			reverse_kmer(adjkmer, adjrevkmer, k/2);
			strcpy(revreference, revkmer);
			strcat(revreference, adjrevkmer);

			fullhash = hash(reference, k);
			adjhash = hash(adjkmer, k/2);
			count += (unsigned long)dbg->nodes[i]->out_kstep[adjhash]->count;
			count_input += (unsigned long)dbg->nodes[i]->out_kstep[adjhash]->input_count;

			revhash = hash(revkmer, k/2);
			revadjhash = hash(adjrevkmer, k/2);
			if(!is_palyndrome(reference, revreference)) {
				count += (unsigned long)dbg->nodes[revhash]->out_kstep[revadjhash]->count;
				count_input += (unsigned long)dbg->nodes[revhash]->out_kstep[revadjhash]->input_count;
				for(j=0; j<k; j++) {
					psm[fullhash][get_base_index(reference[j])][j] += dbg->nodes[revhash]->out_kstep[revadjhash]->count;
				}
			}

			for(j=0; j<k; j++) {
				psm[fullhash][get_base_index(reference[j])][j] += dbg->nodes[i]->out_kstep[adjhash]->count;
			}
			substitute_all(adjkmer, substituted_adj, k/2); //All substitution of second half

			//Iterate over substitution
			for(j=0; j<expected_sub; j++) {
				substituted_node = dbg->nodes[hash(substituted[j], k/2)];
				reverse_kmer(substituted_node->seq, revkmer, k/2);
				rev_substituted_node = dbg->nodes[hash(revkmer, k/2)];

				//Iterate over adjacent k-mer
				count += (unsigned long)substituted_node->out_kstep[adjhash]->count;
				count_input += (unsigned long)substituted_node->out_kstep[adjhash]->input_count;
				strcpy(support, substituted[j]);
				strcat(support, adjkmer);
				strcpy(revsupport, revkmer);
				strcat(revsupport, adjrevkmer);
				if(!is_palyndrome(support, revsupport)) {
					count += (unsigned long)rev_substituted_node->out_kstep[revadjhash]->count;
					count_input += (unsigned long)rev_substituted_node->out_kstep[revadjhash]->input_count;
					for(h=0; h<k; h++) {
						psm[fullhash][get_base_index(support[h])][h] += rev_substituted_node->out_kstep[revadjhash]->count;
					}
				}

				for(h=0; h<k; h++) {
					psm[fullhash][get_base_index(support[h])][h] += substituted_node->out_kstep[adjhash]->count;
				}

				for(h=0; h<expected_sub; h++) {
					subhash = hash(substituted_adj[h], k/2);
					count += (unsigned long)substituted_node->out_kstep[subhash]->count;
					count_input += (unsigned long)substituted_node->out_kstep[subhash]->input_count;

					reverse_kmer(dbg->nodes[subhash]->seq, adjrevkmer, k/2);
					revsubhash = hash(revkmer, k/2);
					strcpy(support, substituted[j]);
					strcat(support, substituted_adj[h]);
					reverse_kmer(substituted[j], revkmer, k/2);
					strcpy(revsupport, revkmer);
					strcat(revsupport, adjrevkmer);
					if(!is_palyndrome(support, revsupport)) {
						count += (unsigned long)rev_substituted_node->out_kstep[revsubhash]->count;
						count_input += (unsigned long)rev_substituted_node->out_kstep[revsubhash]->input_count;
						for(x=0; x<k; x++) {
							psm[fullhash][get_base_index(support[x])][x] += rev_substituted_node->out_kstep[revsubhash]->count;
						}
					}

					for(x=0; x<k; x++) {
						psm[fullhash][get_base_index(support[x])][x] += substituted_node->out_kstep[subhash]->count;
					}
				}
			}
			//Missing counts: first half fixed, second half substituted
			reverse_kmer(kmer, revkmer, k/2);
			for(h=0; h<expected_sub; h++) {
				subhash = hash(substituted_adj[h], k/2);
				count += (unsigned long)dbg->nodes[i]->out_kstep[subhash]->count;
				count_input += (unsigned long)dbg->nodes[i]->out_kstep[subhash]->input_count;

				reverse_kmer(dbg->nodes[subhash]->seq, adjrevkmer, k/2);
				revsubhash = hash(adjrevkmer, k/2);
				strcpy(support, kmer);
				strcat(support, substituted_adj[h]);
				strcpy(revsupport, revkmer);
				strcat(revsupport, adjrevkmer);
				if(!is_palyndrome(support, revsupport)) {
					count += (unsigned long)dbg->nodes[revhash]->out_kstep[revsubhash]->count;
					count_input += (unsigned long)dbg->nodes[revhash]->out_kstep[revsubhash]->input_count;
					for(x=0; x<k; x++) {
						psm[fullhash][get_base_index(support[x])][x] += dbg->nodes[revhash]->out_kstep[revsubhash]->count;
					}
				}

				for(x=0; x<k; x++) {
					psm[fullhash][get_base_index(support[x])][x] += dbg->nodes[i]->out_kstep[subhash]->count;
				}
			}
			counts[hash(reference, k)] = count;
			counts_input[hash(reference, k)] = count_input;
			total += count;
			total_input += count_input;
		}
	}
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Counting ended\n");
	fprintf(stdout, "                  Counted %lu total k-mer for IP\n", total);
	fprintf(stdout, "                  Counted %lu total k-mer for Input\n", total_input);

	//// OUTPUT
	if(c) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: count\n");
		strncpy(out_file, input_file, BUFFER);
		sprintf(buf, ".k%d", k);
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
		for(i=0; i<nodes*nodes; i++) {
			rev_hash(i, k, reference);
			freq = (double)counts[i]/(double)total;
			freq_input = (double)counts_input[i]/(double)total_input;
			diff = freq/freq_input;
			diff_log2 = log2(diff);
			fprintf(fp, "%s\t%lu\t%lf\t%lu\t%lf\t%lf\t%lf\n", reference, counts[i], freq, counts_input[i], freq_input, diff, diff_log2);
		}
		fclose(fp);
	}

	if(psm_arg) {
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: position specific matrix\n");
		FILE * fp_psm;
		strncpy(out_file, input_file, BUFFER);
		sprintf(buf, ".k%d", k);
		strcat(out_file, buf);
		strcat(out_file, ".psm");
		if( !(fp_psm = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		for(i=0; i<nodes*nodes; i++) {
			rev_hash(i, k, reference);
			fprintf(fp_psm, "%s\n", reference);
			fprintf(fp_psm, "A");
			for(j=0; j<k; j++) {
				fprintf(fp_psm, "\t%d", psm[i][0][j]);
			}
			fprintf(fp_psm, "\n");
			fprintf(fp_psm, "C");
			for(j=0; j<k; j++) {
				fprintf(fp_psm, "\t%d", psm[i][1][j]);
			}
			fprintf(fp_psm, "\n");
			fprintf(fp_psm, "G");
			for(j=0; j<k; j++) {
				fprintf(fp_psm, "\t%d", psm[i][2][j]);
			}
			fprintf(fp_psm, "\n");
			fprintf(fp_psm, "T");
			for(j=0; j<k; j++) {
				fprintf(fp_psm, "\t%d", psm[i][3][j]);
			}
			fprintf(fp_psm, "\n");
		}
		fclose(fp_psm);
	}

	if(g) {
		char path[BUFFER+1];
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
		fprintf(stdout, "Generating output: graph\n");
		strcpy(path, input_file);
		sprintf(buf, ".k%d.graph", k);
		strcat(path, buf);
		strcpy(out_file, path);
		strcat(out_file, ".nodes");
		FILE * fp_nodes;
		if( !(fp_nodes = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		strcpy(out_file, path);
		strcat(out_file, ".edges.out");
		FILE * fp_edges_out;
		if( !(fp_edges_out = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		strcpy(out_file, path);
		strcat(out_file, ".edges.out");
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











//////////////////////// FUNCTIONS

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
