#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"

#ifdef SILENT
	#define printf(...)
#endif


/////////////// DEFINES
#define BUFFER 256 //Input buffer
#define ARGS 2
#define IN_FORMAT 1
#define K_ARG 2
#define L_ARG 4
#define IN_FILE 6
#define ARGS 2 //min args
#define ARGS_BUF 50

#define FASTA 1
#define FASTQ 2


/////////////// DATA STRUCTRUES

typedef struct graph_s {
	node_t ** nodes; //hash table of nodes
	edge_t ** edges; //hash table of edges
} graph_t;

/////////////// GLOBAL
int k = 6;
int l = 34;
char ** kmers;
char * subkmers[2];

/////////////// PROTOTYPES

int de_bruijn_ize(char *, graph_t *);
void extract_kmers(char *);
void extract_subkmers();

int contains(char *, char);
int hash(char *);

/////////////// MAIN
int main (int argc, char * argv[]) {
	//// Check args
	int input_format = 0;
	int l_arg = L_ARG;
	int in_file = IN_FILE;
	char input_file[ARGS_BUF+7]; //plus extension
	char nodes_file[ARGS_BUF+7];
	char edges_file[ARGS_BUF+7];
	argc -= 1;
	if (argc < ARGS) {
		fprintf(stdout, "Usage: (--fasta|--fastq) [-k K (default 6)] [-l reads_len (default 34)] input_file_no_ext\n");
		return 1;
	} else {
		if ( strcmp(argv[K_ARG], "-k") == 0 ) {
			k = atoi(argv[K_ARG + 1]);
		} else {
			in_file -= 2;
			l_arg -= 2;
		}
		if ( strcmp(argv[l_arg], "-l") == 0 ){
			l = atoi(argv[l_arg + 1]);
		} else {
			in_file -= 2;
		}
		strncpy(input_file, argv[in_file], ARGS_BUF);
		strcpy(nodes_file, input_file);
		strcat(nodes_file, ".nodes");
		strcpy(edges_file, input_file);
		strcat(edges_file, ".edges");
		if (strcmp("--fasta", argv[IN_FORMAT]) == 0) {
			input_format = FASTA;
			strcat(input_file, ".fa");
		} else if (strcmp("--fastq", argv[IN_FORMAT]) ==  0) {
			input_format = FASTQ;
			strcat(input_file, ".fastq");
		} else {
			fprintf(stdout, "Usage: (--fasta|--fastq) [-k K (default 34)] input_file_no_ext\n");
			return 1;
		}
	}

	//// Variables

	char buf[BUFFER+1];
	FILE * in_fp;
	int i;

	char * read;
	if ( !(read = (char*)malloc(sizeof(char) * (l+1))) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}

	if ( !(kmers = (char**)malloc(sizeof(char*) * (l-k+1) )) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}
	for (i=0; i<(l-k+1); i++) {
		if ( !(kmers[i] = (char*)malloc(sizeof(char) * (k+1) )) ) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return 1;
		}
	}

	for(i=0; i<2; i++) {
		if( !(subkmers[i] = (char*)malloc(sizeof(char) * k)) ) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return 1;
		}
	}

	double nodes = pow((double)4, (double)(k-1));
	double edges = pow((double)4, (double)k);

	graph_t dbg;
	if ( !(dbg.nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}
	if ( !(dbg.edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}
	for (i=0; i<nodes; i++) {
		dbg.nodes[i] = NULL;
		dbg.edges[i] = NULL;
	}
	for( ; i<edges; i++) {
		dbg.edges[i] = NULL;
	}


	//// Opening file
	fprintf(stdout, "Opening: %s\n", input_file);
	if (!(in_fp = fopen(input_file, "r"))) {
		fprintf(stdout, "ERROR: can't open %s\n", input_file);
		return 1;
	}

	//// Reading file line by line
	int skip_line;
	if (input_format == FASTA)
		skip_line = 2;
	else //is FASTQ
		skip_line = 4;

	i=0;
	fgets(buf, BUFFER, in_fp);
	while(!feof(in_fp)) {
		i++;
		if (i==2) {
			strncpy(read, buf, l);
			read[l] = '\0';
			printf("%s\n", read);
			if(de_bruijn_ize(read, &dbg)) {
				return 1;
			}
		}
		if (i==skip_line) {
			i=0;
		}

		fgets(buf, BUFFER, in_fp);
	}

	fclose(in_fp);

	fprintf(stdout, "Generating output files for %s\n", input_file);
	//// Graph computed, proceed to output
	FILE * nodes_fp;
	FILE * edges_fp;

	printf("Opening: %s\n", nodes_file);
	if (!(nodes_fp = fopen(nodes_file, "w+"))) {
		fprintf(stdout, "ERROR: can't open nodes output file\n");
		return 1;
	}
	printf("Opening: %s\n", edges_file);
	if (!(edges_fp = fopen(edges_file, "w+"))) {
		fprintf(stdout, "ERROR: can't open edges output file\n");
		return 1;
	}
	//Printing headers
	fprintf(nodes_fp, "id\tseq\t\n");
	fprintf(edges_fp, "id\tcount\tfrom_node_id\tto_node_id\n");
	for(i=0; i<nodes; i++) {
		if( dbg.nodes[i] ) {
			fprintf(nodes_fp, "%d\t%s\n", (dbg.nodes[i])->id, (dbg.nodes[i])->seq);
		}
		if( dbg.edges[i] ) {
			fprintf(edges_fp, "%d\t%d\t%d\t%d\n",  (dbg.edges[i])->id, (dbg.edges[i])->count, ((dbg.edges[i])->from)->id, ((dbg.edges[i])->to)->id);
		}
	}
	fclose(nodes_fp);
	printf("Remaining edges\n");
	for(i=nodes ; i<edges; i++) {
		if( dbg.edges[i] ) {
			fprintf(edges_fp, "%d\t%d\t%d\t%d\n",  (dbg.edges[i])->id, (dbg.edges[i])->count, ((dbg.edges[i])->from)->id, ((dbg.edges[i])->to)->id);
		}
	}
	fclose(edges_fp);


	return 0;
}



/////////////// FUNCTIONS

int de_bruijn_ize(char * seq, graph_t * graph) {

	extract_kmers(seq);

	int i;

	int hashed[2];
	int edge_id;

	node_t * n0;
	node_t * n1;
	edge_t * e;

	for (i=0; i<l-k+1; i++) {
		//For each kmer
		if (!contains(kmers[i], 'N')) {
			extract_subkmers(kmers[i]);
			printf("\t\t%s\t%s\n", subkmers[0], subkmers[1]);
			hashed[0] = hash(subkmers[0]);
			if( !(n0 = graph->nodes[hashed[0]]) ) {
				if( !(n0 = create_node(hashed[0], subkmers[0])) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				if ( !(graph->nodes[hashed[0]] = (node_t*)malloc(sizeof(node_t))) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				graph->nodes[hashed[0]] = n0;
				printf("added(n0) - ");
			}

			hashed[1] = hash(subkmers[1]);
			if( !(n1 = graph->nodes[hashed[1]]) ) {
				if( !(n1 = create_node(hashed[1], subkmers[1])) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				if ( !(graph->nodes[hashed[1]] = (node_t*)malloc(sizeof(node_t))) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				graph->nodes[hashed[1]] = n1;
				printf("added(n1) - ");
			}
			printf("n0: %d - %s\tn1: %d - %s\n", n0->id, n0->seq, n1->id, n1->seq);
			edge_id = hash(kmers[i]);

			if (!( e = graph->edges[edge_id] )) {
				//Create edge
				if( !(e = create_edge(n0, n1, edge_id) ) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				if ( !(graph->edges[edge_id] = (edge_t*)malloc(sizeof(edge_t))) ) {
					fprintf(stdout, "ERROR: couldn't allocate memory\n");
					return 1;
				}
				graph->edges[edge_id] = e;
				n0 = add_out_edges(n0, e);
				if (!n0) {
					fprintf(stdout, "ERROR: couldn't allocate\n");
					return 1;
				}
				n1 = add_in_edges(n1, e);
				if (!n1) {
					fprintf(stdout, "ERROR: couldn't allocate\n");
					return 1;
				}
				printf("created - ");
			} else {
				//Increment count
				e->count = (e->count) + 1;
			}
			printf("edge: %d\tcount: %d\n", e->id, e->count);
		}

	}
	return 0;
}

void extract_kmers(char * seq) {
	int i, j;
	printf("Extracting kmers from: %s\n", seq);
	for (i=0; i<l-k+1; i++) {
		for(j=0; j<k; j++) {
			kmers[i][j] = seq[i+j];
		}
		kmers[i][j] = '\0';
		printf("%s\n", kmers[i]);
	}
}

void extract_subkmers(char * kmer) {
	int i;
	subkmers[0][0] = kmer[0];
	for(i=1; i<k; i++) {
		subkmers[0][i] = kmer[i];
		subkmers[1][i-1] = kmer[i];
	}
	subkmers[1][i-1] = kmer[i];

	subkmers[0][k-1] = '\0';
	subkmers[1][k-1] = '\0';
}


int contains(char * seq, char c) {
	unsigned int i;
	for(i=0; i<strlen(seq); i++) {
		if(seq[i] == c) {
			return 1;
		}
	}
	return 0;
}

int hash(char * kmer) {
	int hashed = 0;
	unsigned int i;
	int v;
	for(i=0; i < strlen(kmer); i++) {
		switch(kmer[i]) {
			case 'A':
				v = 0;
				break;
			case 'C':
				v = 1;
				break;
			case 'G':
				v = 2;
				break;
			case 'T':
				v = 3;
				break;
			default:
				v = 0;
		}
		hashed = hashed << 2;
		hashed += v;
	}

	return hashed;
}
