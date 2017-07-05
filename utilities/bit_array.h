#ifndef BIT_ARRAY_H
	#define BIT_ARRAY_H

	/* Bit array struct
	 * @field size: bit size
	 * @field real_length: actual length of the array of char
	 * @fiel bits: pointer to array of char (used as bit array)
	 *
	 */
	typedef struct bit_array_s {
		int size;
		int real_length;
		char * bits;
	} bit_array_t;

	/* Init an all zeros bit array
	 * @param size: bit size
	 * @return b: bit_array struct (or NULL if error)
	 *
	 */
	bit_array_t * init_bit_array(int size);

	/* Clear a bit array
	 * @param b: array to free
	 *
	 */
	void free_bit_array(bit_array_t * b);

	/* Get the i-th bit
	 * @param b: bit array
	 * @param i: index to check
	 * @return a positive value if is 1, 0 otherwise, -1 if error
	 *
	 */
	int get_bit(bit_array_t * b, int i);

	/* Set the i-th bit to 1
	 * @param b: bit array
	 * @param i: index to set
	 * @return b: bit_array struct (or NULL if error)
	 *
	 */
	bit_array_t * set_on_bit(bit_array_t * b, int i);

	/* Set the i-th bit to 0
	 * @param b: bit array
	 * @param i: index to set
	 * @return b: bit_array struct (or NULL if error)
	 *
	 */
	bit_array_t * set_off_bit(bit_array_t * b, int i);

	/* Set all array to zeros
	 * @param b: bit array
	 * @return bit array
	 *
	 */
	bit_array_t * clear_bit_array(bit_array_t * b);

#endif
