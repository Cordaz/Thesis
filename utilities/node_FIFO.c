#include "node_FIFO.h"
#include <stdlib.h>

node_FIFO_t * init_queue(int l) {
	node_FIFO_t * q;
	elem_t * e;
	elem_t * new;
	int i;

	if( !(q = (node_FIFO_t*)malloc(sizeof(node_FIFO_t))) ) {
		return NULL;
	}

	//Init first elem
	if ( !(new = (elem_t*)malloc(sizeof(elem_t))) ) {
		return NULL;
	}
	q->first = new;
	q->last = new;
	//Init all others
	for (i=1; i<l; i++) {
		e = new;
		if ( !(new = (elem_t*)malloc(sizeof(elem_t))) ) {
			return NULL;
		}
		e->next = new;
	}
	new->next = q->last;
	q->last = new;

	return q;
}


node_t * dequeue(node_FIFO_t * q) {
	node_t * n = q->first->n;
	q->first = (q->first)->next;

	return n;
}


node_FIFO_t * enqueue(node_FIFO_t * q, node_t * n) {
	q->last->next->n = n;
	q->last = q->last->next;
	return q;
}
