#ifndef SET_H
	#define SET_H

	typedef struct set_s {
		char ** data;
		int first;
		int last;
		int valid_items;
		int dim;
	} set_t;

	set_t * initialize_set(int, int);

	int is_empty(set_t * );

	int put(set_t *, char *);

	int get(set_t *, char *);

	void free_set(set_t *);

	int is_in(set_t *, char *);

	void clear(set_t *);

#endif
