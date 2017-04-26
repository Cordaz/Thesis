#include <stdlib.h>
#include <string.h>
#include "list.h"

list_t * add(list_t * head, char * str) {
	list_t * t;
	list_t * new;
	if(!head) {
		if( !(head = (list_t*)malloc(sizeof(list_t))) ) {
			return NULL;
		}
		if( !(head->str = (char*)malloc(sizeof(char) * strlen(str))) ) {
			return NULL;
		}
		strcpy(head->str, str);
		head->next = NULL;

		return head;
	}

	t = head;
	while(t->next)
		t = t->next;

	if( !(new = (list_t*)malloc(sizeof(list_t))) ) {
		return NULL;
	}
	if( !(new->str = (char*)malloc(sizeof(char) * strlen(str))) ) {
		return NULL;
	}
	new->next = NULL;
	strcpy(new->str, str);
	t->next = new;

	return head;
}

int search(list_t * head, char * str) {
	list_t * t;
	t = head;
	while(t) {
		if( strcmp(t->str, str) == 0 )
			return 1;
		t = t->next;
	}

	return 0;
}

void free_list(list_t * head) {
	if(!head)
		return;

	free_list(head->next);
	free(head);
}
