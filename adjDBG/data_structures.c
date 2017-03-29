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
	new->in_kstep = NULL;
	new->out_kstep = NULL;

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
	tree_edge_t * new;
	tree_edge_t * cur;
	tree_edge_t * prev;

	if( !(new = (tree_edge_t*)malloc(sizeof(tree_edge_t))) ) {
		return NULL;
	}

	new->e = e;

	prev = NULL;
	cur = n->in_kstep;

	while(cur) {
		prev = cur;
		if (new->e->from->id < cur->e->from->id) {
			cur = cur->left;
		} else {
			cur = cur->right;
		}
	}
	new->p = prev;
	if(prev == NULL) {
		//Was Empty
		n->in_kstep = new;
	} else if (new->e->from->id < prev->e->from->id) {
		prev->left = new;
	} else {
		prev->right = new;
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
	tree_edge_t * new;
	tree_edge_t * cur;
	tree_edge_t * prev;

	if( !(new = (tree_edge_t*)malloc(sizeof(tree_edge_t))) ) {
		return NULL;
	}

	new->e = e;

	prev = NULL;
	cur = n->out_kstep;

	while(cur) {
		prev = cur;
		if (new->e->to->id < cur->e->to->id) {
			cur = cur->left;
		} else {
			cur = cur->right;
		}
	}
	new->p = prev;
	if(prev == NULL) {
		//Was Empty
		n->out_kstep = new;
	} else if (new->e->to->id < prev->e->to->id) {
		prev->left = new;
	} else {
		prev->right = new;
	}

	return n;

}

tree_edge_t * search_tree(tree_edge_t * t, int key) {
	if (t == NULL) {
		return NULL;
	}
	if (t->e->to->id == key) {
		return t;
	}
	if (key < t->e->to->id) {
		return search_tree(t->left, key);
	}
	return search_tree(t->right, key);
}

edge_t * exist_edge(node_t * n0, node_t * n1) {
	tree_edge_t * t = n0->out_kstep;
	t = search_tree(t, n1->id);

	if (t)
		return t->e;
	return NULL;
}


edge_t * create_edge(node_t * from, node_t * to, int hashed) {
	edge_t * new;
	if ( ( new = (edge_t *)malloc(sizeof(edge_t)) ) ) {
		new->id = hashed;
		new->count = 0;
		new->input_count = 0;
		new->from = from;
		new->to = to;
	}

	return new;
}
