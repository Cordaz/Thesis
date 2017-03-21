#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"


//// FUNCTIONS
node_t * create_node(int hashed, char * seq) {
	node_t * new;

	if( (new = (node_t *) malloc(sizeof(node_t)) ) ) {
		new->id = hashed;
		strcpy(new->seq, seq);
		new->in = NULL;
		new->out = NULL;
	}

	return new;
}


node_t * add_in_edges(node_t * n, edge_t * e) {
	list_edge_t * cur, * prev;
	cur = n->in;
	prev = NULL;
	while (cur) {
		prev = cur;
		cur = (n->in)->next;
	}

	if(cur) {
		//This should be NULL if here, ERROR
		fprintf(stdout, "ERROR: NULL pointer expected\n");
		exit(1);
	}

	list_edge_t * new;

	//Build new list entry
	if( !( new = (list_edge_t *)malloc(sizeof(list_edge_t)) ) ) {
		fprintf(stdout, "ERROR: couldn't allocate\n");
		exit(1);
	}
	new->e = e;
	new->next = NULL;

	if( !prev ) {
		//Non existing list, add reference in node
		n->in = new;
	} else {
		prev->next = new;
	}

	return n;

}


node_t * add_out_edges(node_t * n, edge_t * e) {
	list_edge_t * cur, * prev;
	cur = n->out;
	prev = NULL;
	while (cur) {
		prev = cur;
		cur = (n->out)->next;
	}

	if(cur) {
		//This should be NULL if here, ERROR
		fprintf(stdout, "ERROR: NULL pointer expected\n");
		exit(1);
	}

	list_edge_t * new;

	//Build new list entry
	if( !( new = (list_edge_t *)malloc(sizeof(list_edge_t)) ) ) {
		fprintf(stdout, "ERROR: couldn't allocate\n");
		exit(1);
	}
	new->e = e;
	new->next = NULL;

	if( !prev ) {
		//Non existing list, add reference in node
		n->out = new;
	} else {
		prev->next = new;
	}

	return n;

}


edge_t * create_edge(node_t * from, node_t * to, int hashed) {
	edge_t * new;
	if ( ( new = (edge_t *)malloc(sizeof(edge_t)) ) ) {
		new->id = hashed;
		new->count = 1;
		new->from = from;
		new->to = to;
	}

	return new;
}
