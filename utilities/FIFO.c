#include "FIFO.h"
#include <stdlib.h>

fifo_t * init_queue(int l) {
	fifo_t * q;
	elem_t * e;
	elem_t * new;
	int i;

	if( !(q = (fifo_t*)malloc(sizeof(fifo_t))) ) {
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


node_t * dequeue(fifo_t * q) {
	node_t * n = q->first->n;
	q->first = (q->first)->next;

	return n;
}


fifo_t * enqueue(fifo_t * q, node_t * n) {
	q->last->next->n = n;
	q->last = q->last->next;
	return q;
}
