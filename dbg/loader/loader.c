#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include "../data_structures.h"

graph_t * alloc_graph(double nodes, double edges) {
	graph_t * dbg;
	int i;

	if ( !(dbg = (graph_t*)malloc(sizeof(graph_t))) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}

	if ( !(dbg->nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}
	if ( !(dbg->edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}
	for (i=0; i<nodes; i++) {
		dbg->nodes[i] = NULL;
		dbg->edges[i] = NULL;
	}
	for( ; i<edges; i++) {
		dbg->edges[i] = NULL;
	}

	return dbg;
}


graph_t * load_nodes(graph_t * dbg, char * path) {
	FILE * fp;
	char buf[64]; //Reading buffer
	char * args; //Args buffer
	node_t * n;

	if ( !(fp = fopen(path, "r")) ) {
		fprintf(stdout, "ERROR: can't open %s\n", path);
		return NULL;
	}

	fgets(buf, 64, fp);
	//Header read, skip
	fgets(buf, 64, fp);
	//Actually read the first line
	int l; //Store length
	while(!feof(fp)) {
		l =  strlen(buf);
		buf[l-1] = '\0'; //remove \n
		if( !(n = (node_t *)malloc(sizeof(node_t))) ) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return NULL;
		}
		//Assume the string is correctly formatted
		if(args = strtok(buf, "\t"))
			n->id = atoi(args);
		if(args = strtok(NULL, "\t"))
			strcpy(n->seq, args);

		dbg->nodes[n->id] = n;

		fgets(buf, 64, fp);
	}

	return dbg;
}


graph_t * load_edges(graph_t * dbg, char * path) {
	FILE * fp;
	char buf[64]; //Reading buffer
	char * args; //Args buffer
	edge_t * e;

	if ( !(fp = fopen(path, "r")) ) {
		fprintf(stdout, "ERROR: can't open %s\n", path);
		return NULL;
	}

	fgets(buf, 64, fp);
	//Header read, skip
	fgets(buf, 64, fp);
	//Actually read the first line
	int l; //Store length
	while(!feof(fp)) {
		l =  strlen(buf);
		buf[l-1] = '\0'; //remove \n
		if( !(e = (edge_t *)malloc(sizeof(edge_t))) ) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return NULL;
		}
		//Assume the string is correctly formatted
		if(args = strtok(buf, "\t"))
			e->id = atoi(args);
		if(args = strtok(NULL, "\t"))
			e->count = atoi(args);
		if(args = strtok(NULL, "\t"))
			e->from = dbg->nodes[atoi(args)];
		if(args = strtok(NULL, "\t"))
			e->to = dbg->nodes[atoi(args)];

		dbg->edges[e->id] = e;

		//Update correlated node
		if (!add_out_edges(e->from, e)) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return NULL;
		}
		if (!add_in_edges(e->to, e)) {
			fprintf(stdout, "ERROR: couldn't allocate memory\n");
			return NULL;
		}

		fgets(buf, 64, fp);
	}

	return dbg;
}
