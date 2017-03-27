#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"


//// FUNCTIONS
node_t * create_node(int hashed, char * seq) {
	node_t * new;

	if( !(new = (node_t *) malloc(sizeof(node_t)) ) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}

	new->id = hashed;
	if ( !(new->seq = (char*)malloc(sizeof(char) * strlen(seq))) ) {
		fprintf(stdout, "ERROR: couldn't allocate memory\n");
		return NULL;
	}
	strcpy(new->seq, seq);
	new->in = NULL;
	new->out = NULL;

	return new;
}


node_t * add_in_edges(node_t * n, edge_t * e) {
	list_edge_t * cur, * prev;
	cur = n->in;
	prev = NULL;
	while (cur) {
		prev = cur;
		cur = cur->next;
	}

	list_edge_t * new;

	//Build new list entry
	new = (list_edge_t *)malloc(sizeof(list_edge_t));
	if( !new ) {
		return NULL;
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

node_t * add_in_kstep_edges(node_t * n, edge_t * e) {
	list_edge_t * cur, * prev;
	cur = n->in_kstep;
	prev = NULL;
	while (cur) {
		prev = cur;
		cur = cur->next;
	}

	list_edge_t * new;

	//Build new list entry
	new = (list_edge_t *)malloc(sizeof(list_edge_t));
	if( !new ) {
		return NULL;
	}
	new->e = e;
	new->next = NULL;

	if( !prev ) {
		//Non existing list, add reference in node
		n->in_kstep = new;
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
		cur = cur->next;
	}

	list_edge_t * new;

	//Build new list entry
	new = (list_edge_t *)malloc(sizeof(list_edge_t));
	if( !new ) {
		return NULL;
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


node_t * add_out_kstep_edges(node_t * n, edge_t * e) {
	list_edge_t * cur, * prev;
	cur = n->out_kstep;
	prev = NULL;
	while (cur) {
		prev = cur;
		cur = cur->next;
	}

	list_edge_t * new;

	//Build new list entry
	new = (list_edge_t *)malloc(sizeof(list_edge_t));
	if( !new ) {
		return NULL;
	}
	new->e = e;
	new->next = NULL;

	if( !prev ) {
		//Non existing list, add reference in node
		n->out_kstep = new;
	} else {
		prev->next = new;
	}

	return n;

}

int update_edge(edge_t * e) {
	edge_t * t1, * t2;
	node_t * from = e->from;
	node_t * to = e->to;

	if (t1 = exist_edge(from, e, "out")) {
		t1->count += 1;
	} else if (t2 = exist_edge(to, e, "in")) {
		t2->count += 1;
	}
	if(!t1) {
		if (!add_out_kstep_edges(from, e)) {
			return 1;
		}
	}

	if (!t2) {
		if(!add_in_kstep_edges(to, e)) {
			return 1;
		}
	}
	return 0;
}

edge_t * exist_edge(node_t * n, edge_t * e, char * type) {
	list_edge_t * t;
	if( strcmp(type, "out") == 0 ) {
		t = n->out_kstep;
		while(t) {
			if (t->e->to->id == e->to->id && t->e->from->id == e->from->id) {
				return t->e;
			}
			t = t->next;
		}
	} else if ( strcmp(type, "in") == 0 ) {
		t = n->in_kstep;
		while(t) {
			if (t->e->to->id == e->to->id && t->e->from->id == e->from->id) {
				return t->e;
			}
			t = t->next;
		}
	}

	return NULL;
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
