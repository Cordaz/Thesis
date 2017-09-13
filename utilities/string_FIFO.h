/*
@author Andrea Corneo
@year 2017

Implementation of a string FIFO queue.
This implementation takes into account the allocation time and is intended to work with limited size in order to reuse
the same memory allocated.

*/

#ifndef STRING_FIFO_H
	#define STRING_FIFO_H

	typedef struct string_FIFO_s {
		char ** data;
		int first;
		int last;
		int valid_items;
		int dim;
	} string_FIFO_t;

	string_FIFO_t * initialize_set(int, int);

	int is_empty(string_FIFO_t * );

	int put(string_FIFO_t *, char *);

	int get(string_FIFO_t *, char *);

	void free_set(string_FIFO_t *);

	int is_in(string_FIFO_t *, char *);

	void clear(string_FIFO_t *);

#endif
