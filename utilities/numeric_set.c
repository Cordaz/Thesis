#include "numeric_set.h"
#include <malloc.h>
#include <string.h>

numeric_set_t * initialize_set(int dim) {
	int i;

	numeric_set_t * set;

	if( !(set = (numeric_set_t*)malloc(sizeof(numeric_set_t))) ) {
		return NULL;
	}

	if( !(set->data = (int*)malloc(sizeof(int) * dim)) ) {
		return NULL;;
	}

	set->valid_items = 0;
	set->dim = dim;
	set->first = 0;
	set->last = 0;

	return set;
}

int is_empty(numeric_set_t * q) {
	return q->valid_items == 0 ? 1 : 0;
}

int putInt(numeric_set_t * q, int x) {
	if(q->valid_items == q->dim)
		return -1; //Full

	q->data[q->last] = x;
	q->valid_items++;
	q->last = (q->last + 1) /*% q->dim*/;

	return 0;
}

int getInt(numeric_set_t * q, int * px) {
	if(q->valid_items == 0)
		return -1; //Empty

	*px = q->data[q->first];
	q->first = (q->first + 1) /*% q->dim*/;
	q->valid_items--;

	return 0;
}

void free_set(numeric_set_t * q) {
	free(q->data);
	free(q);
}

int is_in(numeric_set_t * q, int x) {
	if(q->valid_items == 0)
		return 0;

	int i;

	for(i=q->first; i<q->last; i++) {
		if( q->data[i] == x ) {
			return 1;
		}
	}

	return 0;
}

void clear(numeric_set_t * q) {
	q->first = 0;
	q->last = 0;
	q->valid_items = 0;
}
