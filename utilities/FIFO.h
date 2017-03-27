#ifndef FIFO_H
	#define FIFO_H

	#include "../dbg/data_structures.h"

	typedef struct elem_s {
		node_t * n;
		struct elem_s * next;
	} elem_t;

	typedef struct fifo_s {
		elem_t * first;
		elem_t * last;
	} fifo_t;

	fifo_t * init_queue(int);

	node_t * dequeue(fifo_t *);

	fifo_t * enqueue(fifo_t *, node_t *);

#endif
