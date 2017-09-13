/*
@author Andrea Corneo
@year 2017

Library for distance functions.
The functions contained may depend on what we needed.

*/

#ifndef DISTANCE_H
	#define DISTANCE_H

	/* Return the Edit distance (up to the specified limit) of the two kmers provided
	 * @param kmer1: first kmer
	 * @param kmer2: second kmer
	 * @param k: length of kmer
	 * @param limit: limit of distance computation
	 * @return d(kmer1, kmer2): an integer between 0 and limit (included), -1 if error
	 *
	 * @require len(kmer1) = len(kmer2) = k (not verified; '\0' excluded
	 *
	 */
	int dist(char * kmer1, char * kmer2, int k, int limit);

#endif
