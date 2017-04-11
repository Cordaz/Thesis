#ifndef NUMERIC_SET_H
	#define NUMERIC_SET_H

	typedef struct numeric_set_s {
		int * data;
		int first;
		int last;
		int valid_items;
		int dim;
	} numeric_set_t;

	numeric_set_t * initialize_set(int);

	int is_empty(numeric_set_t * );

	int putInt(numeric_set_t *, int);

	int getInt(numeric_set_t *, int *);

	void free_set(numeric_set_t *);

	int is_in(numeric_set_t *, int);

	void clear(numeric_set_t *);

#endif
