/**
 * @file bit_reverse.hpp
 * @author Hugues Almorin (hugues.almorin@arelis.com)
 * @brief Bit reverse common tool for FFT
 * @version 0.0.0
 * @date 2021-06-14
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef BIT_REVERSE_HPP_
#define BIT_REVERSE_HPP_

#include <cstdint>
#include <cmath>

/**
 * @brief Reverses bits of input
 * 
 * @tparam M input bit length
 * @param input Input value
 * @return unsigned int Bit reversed output
 */
template<int M>
unsigned int reverse_bits(unsigned int input) {
	int i, rev = 0;
	rev_bits_loop : for (i = 0; i < M; i++) {
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

/**
 * @brief Bit reverse permutation on a vector. In-place
 * 
 * @tparam DTYPE Data type
 * @tparam SIZE Vector size
 * @tparam M Vector bit length
 * @param X_R Vector real
 * @param X_I Vector imag
 */
template<class DTYPE = float, int32_t SIZE, int8_t M>
void bit_reverse(DTYPE X_R[], DTYPE X_I[]) {
	unsigned int reversed;
	unsigned int i;
	DTYPE temp;

	br_loop : for (i = 0; i < SIZE; i++) {
		reversed = reverse_bits<M>(i); // Find the bit reversed index
		if (i <= reversed) {
			// std::cout << "Swapping " << i << " and " << reversed << " sample" << std::endl;
			// Swap the real values
			temp          = X_R[i];
			X_R[i]        = X_R[reversed];
			X_R[reversed] = temp;

			// Swap the imaginary values
			temp          = X_I[i];
			X_I[i]        = X_I[reversed];
			X_I[reversed] = temp;
		}
	}
}

/**
 * @brief Bit reverse permutation on a vector. Out-of-place
 * 
 * @tparam DTYPE Data type
 * @tparam SIZE Vector size
 * @tparam M Vector bit length
 * @param in_R Input Vector real
 * @param in_I Input Vector imag
 * @param out_R Output Vector real
 * @param out_I Output Vector imag
 */
template<class DTYPE = float, int32_t SIZE, int8_t M>
void bit_reverse(DTYPE in_R[], DTYPE in_I[], DTYPE out_R[], DTYPE out_I[]) {
	// #pragma HLS ARRAY_PARTITION variable=in_R dim=1 factor=2 cyclic
	// #pragma HLS ARRAY_PARTITION variable=in_I dim=1 factor=2 cyclic
	// #pragma HLS ARRAY_PARTITION variable=out_R dim=1 factor=2 block
	// #pragma HLS ARRAY_PARTITION variable=out_I dim=1 factor=2 block

	unsigned int reversed;
	unsigned int i;
	// static char rev_addr_lut[64] = {
	// 	0,		32,		16,		48,		 8,		40,		24,		56,
	// 	4,		36,		20,		52,		12,		44,		28,		60,
	// 	2,		34,		18,		50,		10,		42,		26,		58,
	// 	6,		38,		22,		54,		14,		46,		30,		62,
	// 	1,		33,		17,		49,		 9,		41,		25,		57,
	// 	5,		37,		21,		53,		13,		45,		29,		61,
	// 	3,		35,		19,		51,		11,		43,		27,		59,
	// 	7,		39,		23,		55,		15,		47,		31,		63,
	// };

	br_loop : for(i = 0; i < SIZE; i++){
		reversed = reverse_bits<M>(i);
		if(i <= reversed) {
			out_R[i] = in_R[reversed];
			out_R[reversed] = in_R[i];

			out_I[i] = in_I[reversed];
			out_I[reversed] = in_I[i];
		}
	}

	// for (i = 0; i < SIZE; i++) {
	// 	// Swap the real values
	// 	out_R[i] = in_R[rev_addr_lut[i]];

	// 	// Swap the imaginary values
	// 	// out_I[i] = in_I[reversed];
	// 	out_I[i] = in_I[rev_addr_lut[i]];
	// }
}

/**
 * @brief Bit reverse permutation on a vector. Out-of-place from a LUT
 * 
 * @tparam DTYPE Data type
 * @tparam SIZE Vector size
 * @tparam M Vector bit length
 * @param br_table permutation LUT
 * @param in_R Input Vector real
 * @param in_I Input Vector imag
 * @param out_R Output Vector real
 * @param out_I Output Vector imag
 */
template<class DTYPE = float, int32_t SIZE, int8_t M>
void bit_reverse(uint16_t br_table[], DTYPE in_R[], DTYPE in_I[], DTYPE out_R[], DTYPE out_I[]){
	unsigned int i;

	br_loop : for (i = 0; i < SIZE; i++) {
		// Swap the real values
		out_R[i] = in_R[br_table[i]];

		// Swap the imaginary values
		// out_I[i] = in_I[reversed];
		out_I[i] = in_I[br_table[i]];
	}
}

/**
 * @brief Bit reverse permutation on a vector. Out-of-place from a LUT
 * 
 * @tparam DTYPE Data type
 * @tparam SIZE Vector size
 * @tparam M Vector bit length
 * @param br_table permutation LUT
 * @param X_R Vector real
 * @param X_I Vector imag
 */
template<class DTYPE = float, int32_t SIZE, int8_t M>
void bit_reverse(uint16_t br_table[], DTYPE X_R[], DTYPE X_I[]){

	unsigned int reversed;
	unsigned int i;
	DTYPE temp;

	br_loop : for (i = 0; i < SIZE; i++) {
		reversed = br_table[i]; // Find the bit reversed index
		if (i <= reversed) {
			// Swap the real values
			temp          = X_R[i];
			X_R[i]        = X_R[reversed];
			X_R[reversed] = temp;

			// Swap the imaginary values
			temp          = X_I[i];
			X_I[i]        = X_I[reversed];
			X_I[reversed] = temp;
		}
	}
}

#endif // BIT_REVERSE_HPP_
