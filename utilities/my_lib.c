#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "numeric_set.h"

const char bases[4] = {'A', 'C', 'G', 'T'};

int hash(char * kmer, int len) {
	int hashed = 0;
	int i;
	int v;
	for(i=0; i < len; i++) {
		switch(kmer[i]) {
			case 'A':
				v = 0;
				break;
			case 'C':
				v = 1;
				break;
			case 'G':
				v = 2;
				break;
			case 'T':
				v = 3;
				break;
			default:
				v = 0;
		}
		hashed = hashed << 2;
		hashed += v;
	}

	return hashed;
}


void rev_hash(int hash, int k, char * str) {
	int i;
	int mask;
	int c;
	mask = 3; //Get last 2 bit

	for (i=0; i<k; i++) {
		c = hash & mask;
		switch(c) {
			case 0:
				c = 'A';
				break;
			case 1:
				c = 'C';
				break;
			case 2:
				c = 'G';
				break;
			case 3:
				c = 'T';
				break;
			default:
				c = 'N';
		}
		str[k-1-i] = c;
		hash = hash >> 2;
	}
	str[i] = '\0';

}

int contains(char * seq, char c) {
	unsigned int i;
	for(i=0; i<strlen(seq); i++) {
		if(seq[i] == c) {
			return 1;
		}
	}
	return 0;
}


void reverse_kmer(char * kmer, char * rev, int k) {
	int i;
	char c;
	for (i=0; i<k; i++) {
		switch(kmer[k-1-i]) {
			case 'A':
				c = 'T';
				break;
			case 'C':
				c = 'G';
				break;
			case 'G':
				c = 'C';
				break;
			case 'T':
				c = 'A';
				break;
			default:
				c = 'N';
		}
		rev[i] = c;
	}
	rev[i] = '\0';
}


int is_palyndrome(char * kmer, char * rev) {
	return strcmp(kmer, rev) == 0 ? 1 : 0;
}


int get_next_substring(char * seq, int index, int min_l, int * len) {
	int i;
	if(*(seq+index) == '\0')
		return -1;
	//Find a possible start
	i=index;
	while(*(seq+i) == 'N') {
		i++;
	}
	//Can start a new substring
	int j = 0;
	while(seq[i+j] != 'N' && seq[i + j] != '\0') {
		j++;
	}
	if (j >= min_l) {
		*len = j;
		return i;
	}
	return get_next_substring(seq, i+j, min_l, len);
}


int count_substitution(char * ref, char * seq, int k, int max) {
	int i;
	int count = 0;

	for(i=0; i<k && count < max+1; i++) {
		if(ref[i] != seq[i]) {
			count++;
		}
	}

	return count;
}


void substitute_all(char * kmer, char ** substituted, int k) {
	int i, j, h;
	for(i=0; i<k; i++) {
		for(h=0; h<3; h++) {
			strcpy(substituted[i*3+h], kmer);
		}
		h=0;
		for(j=0; j<4; j++) {
			if(kmer[k-i-1] != bases[j]) {
				substituted[i*3+h][k-i-1] = bases[j];
				h++;
			}
		}
	}
}


numeric_set_t * substitute_one_hash(numeric_set_t * set, int hash, int k, int * positional_masks) {
	int i, j;
	int t;

	for(i=0; i<k; i++) {
		t = hash & positional_masks[i];
		for(j=0; j<4; j++) {
			putInt(set, t + (j << 2*i) );
		}
	}

	return set;
}

int * init_positional_masks(int k) {
	int * ms;
	if( !(ms = (int*)malloc(sizeof(int) * k)) ) {
		return NULL;
	}
	int full_ms = (int)pow((double)4, k) - 1;
	int i;
	for(i=0; i<k; i++) {
		ms[i] = full_ms - (3 << 2*i);
	}

	return ms;
}



int get_base_index(char b) {
	switch(b) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			return -1;
	}
}
