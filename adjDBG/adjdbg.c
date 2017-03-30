#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "data_structures.h"
#include "../utilities/my_lib.h"
#include "../utilities/FIFO.h"

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

///////// FUNCTIONS PROTOTYPES

graph_t * build_graph(double, double, int);
int map_read(char *, int, int, graph_t *, fifo_t *);
int map_input_read(char *, int, int, graph_t *, fifo_t *);
node_t * get_successor(node_t *, int, char);

int created_edges;


int main (int argc, char * argv[]) {
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

	int k = K;
	int l = L;
	int o = 0;
	//// BEGIN - PARSING ARGS
	int input_format = 0;
	int l_arg = L_ARG;
	int o_arg = O_ARG;
	int c_file = C_FILE;
	int in_file = IN_FILE;
	char input_file[BUFFER+1];
	char control_file[BUFFER+1];
	char out_file[BUFFER+1];
	argc -= 1;
	if (argc < MIN_ARGS) {
		fprintf(stdout, "Usage: (--fasta|--fastq) [-k K] [-l L] [-o] control_file_no_ext input_file_no_ext\n\n");
		fprintf(stdout, "\t--fasta | --fastq file format\n");
		fprintf(stdout, "\t-k K: length of k-mer (graph built on k/2-mer), default %d\n", K);
		fprintf(stdout, "\t-l L: length of the reads in the input file, default %d\n", L);
		fprintf(stdout, "\t-o outputs the graph in a file named input_file_no_ext.graph\n");
		return 1;
	} else {
		if ( strcmp(argv[K_ARG], "-k") == 0 ) {
			k = atoi(argv[K_ARG + 1]);
		} else {
			in_file -= 2;
			l_arg -= 2;
			o_arg -= 2;
			c_file -= 2;
		}
		if ( strcmp(argv[l_arg], "-l") == 0 ){
			l = atoi(argv[l_arg + 1]);
		} else {
			in_file -= 2;
			o_arg -= 2;
			c_file -= 2;
		}
		if ( strcmp(argv[o_arg], "-o") == 0 ){
			o = 1;
		} else {
			in_file -= 1;
			c_file -= 1;
		}
		strncpy(control_file, argv[c_file], BUFFER);
		strncpy(input_file, argv[in_file], BUFFER);
		strncpy(out_file, input_file, BUFFER);
		strcat(out_file, ".graph");
		if (strcmp("--fasta", argv[IN_FORMAT]) == 0) {
			input_format = FASTA;
			strcat(input_file, ".fa");
			strcat(control_file, ".fa");
		} else if (strcmp("--fastq", argv[IN_FORMAT]) ==  0) {
			input_format = FASTQ;
			strcat(input_file, ".fastq");
			strcat(control_file, ".fastq");
		} else {
			fprintf(stdout, "Usage: (--fasta|--fastq) [-k K (default 34)] [-o] control_file_no_ext input_file_no_ext\n");
			return 1;
		}
	}
	//// END - PARSING ARGS

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
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
	fprintf(stdout, "Empty De Bruijn graph built\n");
	fprintf(stdout, "\t\tcreated %d nodes\n", (int)nodes);
	fprintf(stdout, "\t\tcreated %d 1-step edges\n", (int)edges);
	fprintf(stdout, "\t\tcreated %d %d-step edges\n", (int)(nodes*nodes), k/2);

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
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
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


	created_edges = 0;

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

	//// OUTPUT STATISTICS
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
	fprintf(stdout, "Processing of ChIP-seq complete\n");

	int old_created_edges = created_edges;
	//// MAPPING CONTROL FILE
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
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

	//// OUTPUT STATISTICS
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
	fprintf(stdout, "Processing of Input complete\n");


	if(o) {
		//// OUTPUT
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
		fprintf(stdout, "Generating output\n");
		if( !(fp = fopen(out_file, "w+")) ) {
			fprintf(stdout, "[ERROR] can't open %s\n", out_file);
			return 1;
		}
		fprintf(fp, "node\tin_nodes\tout_nodes\tin_nodes_kstep(node:count-input_count)\tout_nodes_kstep(node:count-input_count)\n");
		for(i=0; i<nodes; i++) {
			fprintf(fp, "%s", dbg->nodes[i]->seq);
			fprintf(fp, "\t%s", dbg->nodes[i]->in[0]->from->seq);
			for(j=1; j<4; j++) {
				fprintf(fp, ";%s", dbg->nodes[i]->in[j]->from->seq);
			}
			fprintf(fp, "\t%s", dbg->nodes[i]->out[0]->to->seq);
			for(j=1; j<4; j++) {
				fprintf(fp, ";%s", dbg->nodes[i]->out[j]->to->seq);
			}
			fprintf(fp, "\t%s:%d-%d", dbg->nodes[i]->in_kstep[0]->from->seq, dbg->nodes[i]->in_kstep[0]->count, dbg->nodes[i]->in_kstep[0]->input_count);
			for(j=1; j<nodes; j++) {
				fprintf(fp, ";%s:%d-%d", dbg->nodes[i]->in_kstep[j]->from->seq, dbg->nodes[i]->in_kstep[j]->count, dbg->nodes[i]->in_kstep[j]->input_count);
			}
			fprintf(fp, "\t%s:%d-%d", dbg->nodes[i]->out_kstep[0]->to->seq, dbg->nodes[i]->out_kstep[0]->count, dbg->nodes[i]->out_kstep[0]->input_count);
			for(j=1; j<nodes; j++) {
				fprintf(fp, ";%s:%d-%d", dbg->nodes[i]->out_kstep[j]->to->seq, dbg->nodes[i]->out_kstep[j]->count, dbg->nodes[i]->out_kstep[j]->input_count);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf(stdout, "[%2d:%2d:%2d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
		fprintf(stdout, "Output generated\n");
	}



	return 0;
}












//////////////////////// FUNCTIONS

int map_read(char * read, int l, int k, graph_t * dbg, fifo_t * q) {
	int i;
	edge_t * e;
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
	edge_t * e;
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
			if( !(e = create_edge(dbg->nodes[i], dbg->nodes[j], i)) ) {
				fprintf(stdout, "[ERROR] couldn't allocate memory\n");
				return NULL;
			}
			dbg->nodes[i]->out_kstep[j] = e;
			dbg->nodes[j]->in_kstep[i] = e;
		}
	}

	return dbg;
}
