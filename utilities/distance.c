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


int dist_shift(char * kmer1, char * kmer2, int k, int limit, int * shift, int max_shift) {
	int i, x;
	int min_d = k; // Maximum distance is limit, just init at higher value (maximum dist among two len-k kmer)
	int d;

	// Positive shift
	for(x=0; x<=max_shift; x++) {
		for(i=0, d=0; i<k-x; i++) {
			if(kmer1[i] != kmer2[i+x]) {
				d++;
			}
			if(d >= limit) { i = k; }
		}
		// Check minimum distance
		if(d < min_d) {
			min_d = d;
			*shift = x;
			// If min_d is now zero return
			if(min_d == 0) { return min_d; }
		}
	}

	// Negative shift
	for(x=1; x<=max_shift; x++) {
		for(i=x, d=0; i<k; i++) {
			if(kmer1[i] != kmer2[i-x]) {
				d++;
			}
			if(d >= limit) { i = k; }
		}
		// Check minimum distance
		if(d < min_d) {
			min_d = d;
			*shift = (-1 * x);
			// If min_d is now zero return
			if(min_d == 0) { return min_d; }
		}
	}

	return min_d;
}

///// TEST
/*
int main(int argc, char * argv[]) {
	int num_test = 9;
	int k = 8;
	int limit = 2;
	int shift;
	int max_shift = 1;
	int d;
	int i;
	//							   Identical / 1 sub	/ 2 sub / shifted 1- / shifted 1+ / shift 1- 1 sub / shifted 1+ 1 sub / shifted 2+ / shifted 2-
	char kmer1[8+1] = "ABCDEFGH";
	char * kmer2[] = { "ABCDEFGH", "ABCDTFGH", "ABCTTFGH", "BCDEFGHT", "TABCDEFG", "BCDETGHT", "TABCTEFG", "TTABCDEF", "CDEFGHTT" };
	// pair dist / shift expected
	int expected_res[9][2] = { {0, 0}, {1, 0}, {2, 0}, {0, -1}, {0, 1}, {1, -1}, {1, 1}, {2, 0}, {2, 0} };

	for(i=0; i<num_test; i++) {
		d = dist_shift(kmer1, kmer2[i], k, limit, &shift, max_shift);
		printf("exp_d=%d\texp_s=%d\ncom_d=%d\tcom_s=%d\npassed=", expected_res[i][0], expected_res[i][1], d, shift);
		if( expected_res[i][0] == d && expected_res[i][1] == shift ) {
			printf("TRUE\n");
		} else {
			printf("FALSE\n");
		}
	}

	return 0;
}
*/
