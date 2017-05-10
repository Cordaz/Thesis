#include "string_FIFO.h"
#include <malloc.h>
#include <string.h>

string_FIFO_t * initialize_set(int dim, int k) {
	int i;

	string_FIFO_t * q;

	if( !(q = (string_FIFO_t*)malloc(sizeof(string_FIFO_t))) ) {
		return NULL;
	}

	if( !(q->data = (char**)malloc(sizeof(char*) * dim)) ) {
		return NULL;;
	}
	for(i=0; i<dim; i++) {
		if( !(q->data[i] = (char*)malloc(sizeof(char) * (k+1))) ) {
			return NULL;
		}
	}

	q->valid_items = 0;
	q->dim = dim;
	q->first = 0;
	q->last = 0;

	return q;
}

int is_empty(string_FIFO_t * q) {
	return q->valid_items == 0 ? 1 : 0;
}

int put(string_FIFO_t * q, char * str) {
	if(q->valid_items == q->dim)
		return -1; //Full

	strcpy(q->data[q->last], str);
	q->valid_items++;
	q->last = (q->last + 1) /*% q->dim*/;

	return 0;
}

int get(string_FIFO_t * q, char * str) {
	if(q->valid_items == 0)
		return -1; //Empty

	strcpy(str, q->data[q->first]);
	q->first = (q->first + 1) /*% q->dim*/;
	q->valid_items--;

	return 0;
}

void free_set(string_FIFO_t * q) {
	int i;
	for(i=0; i<q->dim; i++) {
		free(q->data[i]);
	}
	free(q->data);
	free(q);
}

int is_in(string_FIFO_t * q, char * str) {
	if(q->valid_items == 0)
		return 0;

	int i;

	for(i=q->first; i<q->last; i++) {
		if(strcmp(q->data[i], str) == 0 ) {
			return 1;
		}
	}

	return 0;
}

void clear(string_FIFO_t * q) {
	q->first = 0;
	q->last = 0;
	q->valid_items = 0;
}
