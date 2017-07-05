#include <stdio.h>
#include <malloc.h>
#include <math.h>

#include "bit_array.h"

bit_array_t * init_bit_array(int size) {
	bit_array_t * b;
	if( !(b = (bit_array_t*)malloc(sizeof(bit_array_t))) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return NULL;
	}
	b->size = size;
	b->real_length = (int)ceil( (double)size/(sizeof(char) * 8) );
	if( !(b->bits=(char*)malloc(sizeof(char) * b->real_length)) ) {
		fprintf(stdout, "[ERROR] cannot allocate\n");
		return NULL;
	}
	int i;
	for(i=0; i<b->real_length; i++) {
		b->bits[i] = (char)0x0;
	}

	return b;
}

void free_bit_array(bit_array_t * b) {
	free(b->bits);
	free(b);
}

int get_bit(bit_array_t * b, int i) {
	int index = i / (sizeof(char) * 8);
	if(index >= b->real_length) {
		fprintf(stdout, "[ERROR] out of bounds bit array\n");
		return -1;
	}
	int shift = i % (sizeof(char) * 8);
	unsigned int flag = 1;

	flag = flag << shift;
	return b->bits[index] & flag;
}

bit_array_t * set_on_bit(bit_array_t * b, int i) {
	int index = i / (sizeof(char) * 8);
	if(index >= b->real_length) {
		fprintf(stdout, "[ERROR] out of bounds bit array\n");
		return NULL;
	}
	int shift = i % (sizeof(char) * 8);
	unsigned int flag = 1;

	flag = flag << shift;
	b->bits[index] = b->bits[index] | flag;

	return b;
}

bit_array_t * set_off_bit(bit_array_t * b, int i) {
	int index = i / (sizeof(char) * 8);
	if(index >= b->real_length) {
		fprintf(stdout, "[ERROR] out of bounds bit array\n");
		return NULL;
	}
	int shift = i % (sizeof(char) * 8);
	unsigned int flag = 1;

	flag = flag << shift;
	b->bits[index] = b->bits[index] & (~flag);

	return b;
}

bit_array_t * clear_bit_array(bit_array_t * b) {
	int i;
	for(i=0; i<b->real_length; i++) {
		b->bits[i] = (char)0x0;
	}
	return b;
}
