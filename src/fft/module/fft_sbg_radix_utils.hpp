/**
 * @file fft_sbg_radix_utils.hpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains the utility functions used in the FFT behavioral model
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef FFT_SBG_RADIX_UTILS_HPP_
#define FFT_SBG_RADIX_UTILS_HPP_

#include "../../common/spiral_utils/spiral_utils.hpp"

#include "ap_int.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wshift-count-overflow"

/**
 * @brief Computes twiddle mults for a SW data
 * 
 * @tparam DTYPE data type
 * @tparam TTYPE twiddle factory type
 * @tparam SIZE fft size
 * @tparam RADIX Butterfly radix
 * @tparam SW streaming width
 * @tparam STAGE FFT current stage
 * @param In_R Input data packet (real)
 * @param In_I Input data packet (imag)
 * @param Out_R Output data packet (real)
 * @param Out_I Output data packet (imag)
 * @param Tw_R Twiddle factor table (real)
 * @param Tw_I Twiddle factor table (imag)
 * @param Q Twiddle factor quantization
 */
template<class DTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GSTART>
void fft_sbg_radix_twiddling(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE Out_R[SW], DTYPE Out_I[SW], TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], int STAGE, int Q){
	#pragma HLS INTERFACE ap_ctrl_none port=return
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	
	const unsigned M = ceillogR(SIZE, RADIX);
	constexpr const unsigned log2N = ceillog2(SIZE);
	constexpr const unsigned log2R = ceillog2(RADIX);
	constexpr const unsigned BFMASK = (1 << (log2N-log2R)) - 1;

	ITYPE temp_R;
	ITYPE temp_I;

	static int tw_counter = (SIZE/SW)-1;
	tw_counter = (tw_counter == ((SIZE/SW)-1)) ? 0 : tw_counter + 1;

	for(int bfi = 0; bfi < SW/RADIX; bfi++){
		int divider = ((M-STAGE) * log2R);
		int shiftedbfi = ((tw_counter << divider) & BFMASK);

		// #ifndef __SYNTHESIS__
		// // std::cout << "LOG2N : " << log2N << "; LOG2R" << log2R << std::endl;
   		// std::cout << "ISTAGE : " << STAGE << std::endl;
        // std::cout << "\tBFI : " << bfi << "; DIVIDER : " << (1 << divider) << "; NBF : " << (((tw_counter*SW/RADIX)+bfi) >> divider) << std::endl;
		// #endif

		Out_R[RADIX*bfi] = In_R[RADIX*bfi]; // First twid always with W0
		Out_I[RADIX*bfi] = In_I[RADIX*bfi];

		for(int twi = 1; twi < RADIX; twi++){
			int twid_index = ((((tw_counter*SW/RADIX)+bfi) >> divider)*twi) << divider;
			// #ifndef __SYNTHESIS__
			// std::cout << "\t\tTWI : " << twi << "; TWADDR : " << twid_index << std::endl;
			// #endif

			// cmult<DTYPE, TTYPE, DTYPE, ITYPE>(In_R[RADIX*bfi + twi], In_I[RADIX*bfi + twi], Tw_R[twid_index], Tw_I[twid_index], Out_R[RADIX*bfi + twi], Out_I[RADIX*bfi + twi], Q);

			temp_R = mult_sub<DTYPE, TTYPE, ITYPE>(In_R[RADIX*bfi + twi],Tw_R[twid_index],In_I[RADIX*bfi + twi],-Tw_I[twid_index],Q);
			temp_I = mult_add<DTYPE, TTYPE, ITYPE>(In_I[RADIX*bfi + twi],Tw_R[twid_index],In_R[RADIX*bfi + twi],-Tw_I[twid_index],Q);
			Out_R[RADIX*bfi + twi] = temp_R;
			Out_I[RADIX*bfi + twi] = temp_I;
		}
	}
}

#pragma GCC diagnostic pop

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW>
void fft_sbg_radix_buf_write(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], content_addr<DTYPE, LOG2N, LOG2SW> in_post_switch[SW]){
	for(int i = 0; i < SW; i++){
		#pragma HLS UNROLL
		buf_R[i][in_post_switch[i].addr] = in_post_switch[i].data_R;
		buf_I[i][in_post_switch[i].addr] = in_post_switch[i].data_I;
	}
}

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW>
void fft_sbg_radix_buf_read(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], DTYPE in_from_buff_pre_switch_R[SW], DTYPE in_from_buff_pre_switch_I[SW], ap_uint<LOG2N-LOG2SW+1> ridx_post_switch[SW]){
	for(int i = 0; i < SW; i++){
		#pragma HLS UNROLL
		in_from_buff_pre_switch_R[i] = buf_R[i][ridx_post_switch[i]];
		in_from_buff_pre_switch_I[i] = buf_I[i][ridx_post_switch[i]];
	}
}

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE>
void fft_sbg_radix_digrev_buf_write(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], bool wt_offset, perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config){
	static ap_uint<LOG2N-LOG2SW> j = (SIZE/SW)-1;
	if(j == (SIZE/SW)-1)
		j = 0;
	else
		j++;

	ap_uint<LOG2N-LOG2SW+1> w_addr[SW];
	buf_write_addr_generation<SW, LOG2N, LOG2SW>(j, wt_offset, dig_rev_config.w_addr_bit_seq, w_addr);

	content_addr<DTYPE, LOG2N, LOG2SW> in_pre_switch[SW], in_post_switch[SW];
	combine_addr_data<DTYPE, SW, LOG2N, LOG2SW>(in_pre_switch, w_addr, In_R, In_I);
	switch_network_write<content_addr<DTYPE, LOG2N, LOG2SW>, DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW>(in_pre_switch, in_post_switch, j, dig_rev_config.init_perm_idx, dig_rev_config.w_switch_connection_idx, dig_rev_config.w_switch_control_bit);

	fft_sbg_radix_buf_write<DTYPE, SIZE, SW, LOG2N, LOG2SW>(buf_R, buf_I, in_post_switch);
}

template<class DTYPE, int SIZE, int SW, int GSTART, int LOG2N, int LOG2SW, int STRIDE_PERM_SWITCH_NUM_STAGE>
void fft_sbg_radix_stride_buf_write(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], bool wt_offset, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config){
	static ap_uint<LOG2N-LOG2SW> j = (SIZE/SW)-1;
	if(j == (SIZE/SW)-1)
		j = 0;
	else
		j++;

	ap_uint<LOG2N-LOG2SW+1> w_addr[SW];
	buf_write_addr_generation<SW, LOG2N, LOG2SW>(j, wt_offset, stride_config.w_addr_bit_seq, w_addr);

	content_addr<DTYPE, LOG2N, LOG2SW> in_pre_switch[SW], in_post_switch[SW];
	combine_addr_data<DTYPE, SW, LOG2N, LOG2SW>(in_pre_switch, w_addr, In_R, In_I);
	switch_network_write<content_addr<DTYPE, LOG2N, LOG2SW>, STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW>(in_pre_switch, in_post_switch, j, stride_config.init_perm_idx, stride_config.w_switch_connection_idx, stride_config.w_switch_control_bit);

	fft_sbg_radix_buf_write<DTYPE, SIZE, SW, LOG2N, LOG2SW>(buf_R, buf_I, in_post_switch);
}

template<class DTYPE, int SIZE, int SW, int GSTART, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void fft_sbg_radix_st_dr_buf_read(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], bool rd_offset, bool first_stage, DTYPE Out_R[SW], DTYPE Out_I[SW], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config){
	static ap_uint<LOG2N-LOG2SW> j = (SIZE/SW)-1;
	if(j == (SIZE/SW)-1)
		j = 0;
	else
		j++;

	ap_uint<LOG2N-LOG2SW+1> r_addr[SW], r_addr_post_switch[SW];
	buf_read_addr_generation<SW, LOG2N, LOG2SW>(j, rd_offset, r_addr);

	DTYPE in_from_buff_pre_switch_R[SW], in_from_buff_pre_switch_I[SW];
	fft_sbg_radix_buf_read<DTYPE, SIZE, SW, LOG2N, LOG2SW>(buf_R, buf_I, in_from_buff_pre_switch_R, in_from_buff_pre_switch_I, r_addr);

	if(first_stage){
		switch_network_read<DTYPE, DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW>(in_from_buff_pre_switch_R, in_from_buff_pre_switch_I, Out_R, Out_I, j, dig_rev_config.r_switch_connection_idx, dig_rev_config.r_switch_control_bit);
	}
	else{
		switch_network_read<DTYPE, STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW>(in_from_buff_pre_switch_R, in_from_buff_pre_switch_I, Out_R, Out_I, j, stride_config.r_switch_connection_idx, stride_config.r_switch_control_bit);
	}
}

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int STRIDE_PERM_SWITCH_NUM_STAGE>
void fft_sbg_radix_final_stride_buf_read(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], bool rd_offset, DTYPE Out_R[SW], DTYPE Out_I[SW], perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config){
	static ap_uint<LOG2N-LOG2SW> j = (SIZE/SW)-1;
	if(j == (SIZE/SW)-1)
		j = 0;
	else
		j++;

	ap_uint<LOG2N-LOG2SW+1> r_addr[SW], r_addr_post_switch[SW];
	buf_read_addr_generation<SW, LOG2N, LOG2SW>(j, rd_offset, r_addr);

	DTYPE in_from_buff_pre_switch_R[SW], in_from_buff_pre_switch_I[SW];
	fft_sbg_radix_buf_read<DTYPE, SIZE, SW, LOG2N, LOG2SW>(buf_R, buf_I, in_from_buff_pre_switch_R, in_from_buff_pre_switch_I, r_addr);

	switch_network_read<DTYPE, STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW>(in_from_buff_pre_switch_R, in_from_buff_pre_switch_I, Out_R, Out_I, j, stride_config.r_switch_connection_idx, stride_config.r_switch_control_bit);
}

#endif //FFT_SBG_RADIX_UTILS_HPP
