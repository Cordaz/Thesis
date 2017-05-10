#ifndef FIFO_H
	#define FIFO_H

	#include "data_structures.h"

	typedef struct elem_s {
		node_t * n;
		struct elem_s * next;
	} elem_t;

	typedef struct node_FIFO_s {
		elem_t * first;
		elem_t * last;
	} node_FIFO_t;

	node_FIFO_t * init_queue(int);

	node_t * dequeue(node_FIFO_t *);

	node_FIFO_t * enqueue(node_FIFO_t *, node_t *);

#endif
