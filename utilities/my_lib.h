#ifndef MYLIB_H
	#define MYLIB_H

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

#endif
