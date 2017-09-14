#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"


//// FUNCTIONS
node_t * create_node(int hashed, char * seq, int n) {
	node_t * new;
	int i;

	if( !(new = (node_t *) malloc(sizeof(node_t)) ) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}

	new->id = hashed;
	if ( !(new->seq = (char*)malloc(sizeof(char) * (strlen(seq)+1) )) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}
	strcpy(new->seq, seq);
	for(i=0; i<4; i++) {
		new->in[i] = NULL;
		new->out[i] = NULL;
	}
	if( !(new->in_kstep = (edge_t**)malloc(sizeof(edge_t*) * n)) ) {
		return NULL;
	}
	if( !(new->out_kstep = (edge_t**)malloc(sizeof(edge_t*) * n)) ) {
		return NULL;
	}
	for(i=0; i<n; i++) {
		new->in_kstep[i] = NULL;
		new->out_kstep[i] = NULL;
	}

	return new;
}


edge_t * create_edge(node_t * from, node_t * to, int hashed) {
	edge_t * new;
	if ( ( new = (edge_t *)malloc(sizeof(edge_t)) ) ) {
		new->id = hashed;
		new->count = 0;
		new->kmer_count = 0;
		new->input_count = 1; // pseudo-count
		new->kmer_input_count = 1; // pseudo-count
		new->from = from;
		new->to = to;
	}

	return new;
}
