#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "distance.h"

int dist(char * kmer1, char * kmer2, int k, int limit) {
	int i;
	int d=0;
	for(i=0; i<k; i++) {
		if (kmer1[i] != kmer2[i]) {
			d++;
		}
		if(d >= limit) {
			i = k;
		}
	}
	return d;
}
