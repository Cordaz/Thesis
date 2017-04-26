#ifndef LIST_H
	#define LIST_H

	typedef struct list_s {
		char * str;
		struct list_s * next;
	} list_t;

	list_t * add(list_t *, char *);

	int search(list_t *, char *);

	void free_list(list_t *);

#endif
