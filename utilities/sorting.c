#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "sorting.h"

int * sort_index_decreasing(double * array, int size) {
	int * sorted_index;
	if( !(sorted_index = (int*)malloc(sizeof(int) * size)) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return NULL;
	}

	int i, j;
	int tmp;

	// Init array of index
	for(i=0; i<size; i++) { sorted_index[i] = i; }

	// Bubblesort
	for (i=0; i<size-1; i++) {
		for (j=size-1; j>i; j--) {
			if(array[sorted_index[j]] > array[sorted_index[j-1]]) {
				tmp = sorted_index[j];
				sorted_index[j] = sorted_index[j-1];
				sorted_index[j-1] = tmp;
			}
		}
	}

	return sorted_index;
}

int get_pos(int * sorted, int start, int size, int val) {
	int i;
	for(i=start; i<size; i++) {
		if (sorted[i] == val) return i;
	}
	return -1;
}
