#ifndef MYLIB_H
	#define MYLIB_H

	#include "set.h"

	/*
	 * Hash the string (of length passed after) passed as argument as:
	 * A -> 0 (00)
	 * C -> 1 (01)
	 * G -> 2 (10)
	 * T -> 3 (11)
	 *
	 */
	int hash(char *, int);

	/*
	 * Reverse hash the integer value to a string of length K (second arg)
	 * 0 (00) -> A
	 * 1 (01) -> C
	 * 2 (10) -> G
	 * 3 (11) -> T
	 *
	 */
	void rev_hash(int, int, char *);

	/*
	 * Returns 1 if the string contains the provided char, 0 otherwise
	 *
	 */
	int contains(char *, char);

	/*
	 * Calculates the reverse complement of the prompted kmer and put the result in
	 * the second argument.
	 * The integer arg specify the dimension.
	 *
	 */
	void reverse_kmer(char *, char *, int);

	/*
	 * Returns 1 if the kmer and the complementary reverse are equal (i.e. it's palyndrome), 0 otherwise
	 *
	 */
	int is_palyndrome(char *, char *);

	/*
	 * Get the starting index of the next substring of minimum length that contains no 'N'.
	 * Return -1 if end of string and can't be returned sufficient long string.
	 *
	 */
	int get_next_substring(char *, int, int, int *);

	/*
	 * Count substitution of the two provided k-mer up to length provided
	 * It is bounded to a max+1
	 *
	 */
	int count_substitution(char * ref, char * seq, int k, int max);

	/*
	 * Compute all possible substitution of one nucleotides in the kmer.
	 * Put the result in substituted.
	 *
	 */
	void substitute_all(char *, char **, int);


	set_t * substitute_one(set_t *, char *, int);

	/*
	numeric_set_t * substitute_one_hash(numeric_set_t *, int, int, int *);

	int * init_positional_masks(int);
	*/


	/*
	 * Returns the number corresponding to the base:
	 * A->0, C->1, G->2, T->3, err -1
	 *
	 */
	int get_base_index(char);

#endif
