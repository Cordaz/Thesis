#include <string.h>

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
		switch(kmer[i]) {
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
				c = '\0'; //Used to eventually return error further in the processing, should not happen
		}
		rev[k-1-i] = c;
	}
}
