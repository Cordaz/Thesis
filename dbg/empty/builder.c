#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include "../data_structures.h"
#include "../../utilities/my_lib.h"

#ifdef SILENT
	#define printf(...)
#endif


/////////////// DEFINES
#define BUFFER 256 //Input buffer
#define K_ARG 1

/////////////// GLOBAL
int k = 6;
int mask_hash;

/////////////// MAIN
int main (int argc, char * argv[]) {

	argc -= 1;
	if (argc) {
		if ( (strcmp(argv[K_ARG], "-k") == 0) && (argc == K_ARG + 1)) {
			k = atoi(argv[K_ARG + 1]);
		} else {
			fprintf(stdout, "Argument %s non recognized\n", argv[K_ARG]);
		}
	}

	char out_file[BUFFER+1];
	char kstr[4];
	snprintf(kstr, 4, "%d", k);
	strcpy(out_file, "empty_");
	strcat(out_file, kstr);
	strcat(out_file, ".dbg");

	//// Variables

	int i;

	double nodes = pow((double)4, (double)(k-1));
	double edges = pow((double)4, (double)k);
	mask_hash = (int)nodes - 1;

	graph_t dbg;
	if ( !(dbg.nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}
	if ( !(dbg.edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return 1;
	}

	char * kmer;
	if ( !(kmer = (char*)malloc(sizeof(char) * (k+1))) ) {
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

	int n0_hash, n1_hash;
	node_t * n0, * n1;
	edge_t * e;

	for(i=0; i<edges; i++) {
		//Each i is a coded kmer
		/*
		to retrieve the sequence:
		rev_hash(i, k, kmer);
		*/

		n0_hash = i >> 2;
		n1_hash = i & mask_hash;
		printf("%x\t%x\t%x\n", i, n0_hash, n1_hash);

		if( !(n0 = dbg.nodes[n0_hash]) ) {
			rev_hash(n0_hash, k-1, kmer);
			if( !(n0 = create_node(n0_hash, kmer)) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			if ( !(dbg.nodes[n0_hash] = (node_t*)malloc(sizeof(node_t))) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			dbg.nodes[n0_hash] = n0;
		}

		if( !(n1 = dbg.nodes[n1_hash]) ) {
			rev_hash(n1_hash, k-1, kmer);
			if( !(n1 = create_node(n1_hash, kmer)) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			if ( !(dbg.nodes[n1_hash] = (node_t*)malloc(sizeof(node_t))) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			dbg.nodes[n1_hash] = n1;
		}

		if (!( e = dbg.edges[i] )) {
			//Create edge
			if( !(e = create_edge(n0, n1, i) ) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			if ( !(dbg.edges[i] = (edge_t*)malloc(sizeof(edge_t))) ) {
				fprintf(stdout, "ERROR: couldn't allocate memory\n");
				return 1;
			}
			dbg.edges[i] = e;
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
		}
	}

	//// Graph computed, proceed to output
	FILE * dbg_fp;

	printf("Opening: %s\n", out_file);
	if (!(dbg_fp = fopen(out_file, "w+"))) {
		fprintf(stdout, "ERROR: can't open nodes output file\n");
		return 1;
	}
	//Printing headers
	fprintf(dbg_fp, "edge_id\tfrom[id]\tfrom[seq]\tto[id]\tto[seq]\n");

	for(i=0 ; i<edges; i++) {
		if( dbg.edges[i] ) {
			fprintf(dbg_fp, "%d\t%d\t%s\t%d\t%s\n",  (dbg.edges[i])->id, ((dbg.edges[i])->from)->id, ((dbg.edges[i])->from)->seq, ((dbg.edges[i])->to)->id, ((dbg.edges[i])->to)->seq);
		}
	}
	fclose(dbg_fp);


	return 0;
}
