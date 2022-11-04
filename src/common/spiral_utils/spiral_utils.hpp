// SPIRAL License
//
// Copyright 2017, Carnegie Mellon University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the SPIRAL project.

/**
 * @file spiral_utils.hpp
 * @author Guanglin Xu (guanglix (at) xxxandrew.cmu.edu (delete xxx)), modified by Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains various utility functions used in this SPIRAL project originally written by Guanglin Xu. They have been modified to allow genericity.
 * @version 0.0.0
 * @date 2022-03-01
 * 
 * @license This source is released under the SPIRAL License.
 * 
 */

#ifndef SPIRAL_UTILS_HPP_
#define SPIRAL_UTILS_HPP_

#include "../mult_add.hpp"
#include "stride_perm_num_stage.h"
#include "stride_perm_config.h"
#include "dig_rev_perm_num_stage.h"
#include "dig_rev_perm_config.h"

#include "ap_int.h"

/**
 * @brief Permutation config struct
 * 
 * @tparam NUM_STAGE permutation stages
 * @tparam SW streaming width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 */
template<int NUM_STAGE, int SW, int LOG2N, int LOG2SW>
struct perm_config{
	int init_perm_idx[SW];

	int w_switch_connection_idx[NUM_STAGE][SW];
	int w_switch_control_bit[NUM_STAGE];
	int w_addr_bit_seq[LOG2N-LOG2SW];

	int r_switch_connection_idx[NUM_STAGE][SW];
	int r_switch_control_bit[NUM_STAGE];
};

/**
 * @brief Packer content-address
 * 
 * @tparam DTYPE Data type
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 */
template<class DTYPE, int LOG2N, int LOG2SW>
struct content_addr{
	DTYPE data_R; /**< Real data part */ 
	DTYPE data_I; /**< Imag data part */ 
	ap_uint<LOG2N-LOG2SW+1> addr; /**< Data address*/ 
};

/**
 * @brief Radix 2 butterfly
 * 
 * @tparam DTYPE Data type
 * @param In_R Input data packet (real)
 * @param In_I Input data packet (imag)
 * @param Out_R Output data packet (real)
 * @param Out_I Output data packet (imag)
 */
template<class DTYPE>
void radix_core_2(DTYPE In_R[2], DTYPE In_I[2], DTYPE Out_R[2], DTYPE Out_I[2]){
	Out_R[0] = In_R[0] + In_R[1];
	Out_I[0] = In_I[0] + In_I[1];
	Out_R[1] = In_R[0] - In_R[1];
	Out_I[1] = In_I[0] - In_I[1];
}

/**
 * @brief Radix R butterfly
 * 
 * @tparam DTYPE Data type
 * @tparam RADIX Butterfly radix
 * @param In_R Input data packet (real)
 * @param In_I Input data packet (imag)
 * @param Out_R Output data packet (real)
 * @param Out_I Output data packet (imag)
 */
template<class DTYPE, int RADIX>
void radix_core(DTYPE In_R[RADIX], DTYPE In_I[RADIX], DTYPE Out_R[RADIX], DTYPE Out_I[RADIX]){
	if(RADIX == 2){
		Out_R[0] = In_R[0] + In_R[1];
		Out_I[0] = In_I[0] + In_I[1];
		Out_R[1] = In_R[0] - In_R[1];
		Out_I[1] = In_I[0] - In_I[1];
	}else if(RADIX == 4){
		DTYPE tmp_R[4], tmp_I[4];
		tmp_R[0] = In_R[0] + In_R[2]; // CMPXADD(tmp[0], x[0], x[2]);
		tmp_I[0] = In_I[0] + In_I[2];
		tmp_R[1] = In_R[0] - In_R[2]; // CMPXSUB(tmp[1], x[0], x[2]);
		tmp_I[1] = In_I[0] - In_I[2];
		tmp_R[2] = In_R[1] + In_R[3]; // CMPXADD(tmp[2], x[1], x[3]);
		tmp_I[2] = In_I[1] + In_I[3];
		tmp_R[3] = In_I[1] - In_I[3]; // CMPX_NI_MUL_X_SUB_Y(tmp[3], x[1], x[3]);
		tmp_I[3] = In_R[3] - In_R[1];

		Out_R[0] = tmp_R[0] + tmp_R[2]; // CMPXADD(y[0], tmp[0], tmp[2]);
		Out_I[0] = tmp_I[0] + tmp_I[2];
		Out_R[2] = tmp_R[0] - tmp_R[2]; // CMPXSUB(y[2], tmp[0], tmp[2]);
		Out_I[2] = tmp_I[0] - tmp_I[2];
		Out_R[1] = tmp_R[1] + tmp_R[3]; // CMPXADD(y[1], tmp[1], tmp[3]);
		Out_I[1] = tmp_I[1] + tmp_I[3];
		Out_R[3] = tmp_R[1] - tmp_R[3]; // CMPXSUB(y[3], tmp[1], tmp[3]);
		Out_I[3] = tmp_I[1] - tmp_I[3];
	}else if(RADIX == 8){
		const int Q = 20;
		const float tw = 0.70710678118654757; // 1/sqrtf(2);
		DTYPE s16_R, s17_R, s18_R, s19_R, s20_R;
		DTYPE s16_I, s17_I, s18_I, s19_I, s20_I;
		DTYPE t66_R, t67_R, t68_R, t69_R, t70_R, t71_R, t72_R, t73_R, t74_R, t75_R, t76_R;
		DTYPE t66_I, t67_I, t68_I, t69_I, t70_I, t71_I, t72_I, t73_I, t74_I, t75_I, t76_I;

		t66_R = In_R[0] + In_R[4]; // CMPXADD(t66, X[0], X[4]);
		t66_I = In_I[0] + In_I[4];
		t67_R = In_R[0] - In_R[4]; // CMPXSUB(t67, X[0], X[4]);
		t67_I = In_I[0] - In_I[4];
		t68_R = In_R[1] + In_R[5]; // CMPXADD(t68, X[1], X[5]);
		t68_I = In_I[1] + In_I[5];
		s16_R = mult_litt_f<DTYPE, DTYPE>(tw, ( (In_R[1] - In_R[5]) + (In_I[1] - In_I[5]) ), Q); // CMPX_TWR_NTWR_MUL_X_SUB_Y(s16, 0.70710678118654757, X[1], X[5]);
		s16_I = mult_litt_f<DTYPE, DTYPE>(tw, ( (In_I[1] - In_I[5]) - (In_R[1] - In_R[5]) ), Q);
		t69_R = In_R[2] + In_R[6]; // CMPXADD(t69, X[2], X[6]);
		t69_I = In_I[2] + In_I[6];
		s17_R = In_I[2] - In_I[6]; // CMPX_NI_MUL_X_SUB_Y(s17, X[2], X[6]);
		s17_I = In_R[6] - In_R[2];
		t70_R = In_R[3] + In_R[7]; // CMPXADD(t70, X[3], X[7]);
		t70_I = In_I[3] + In_I[7];
		s18_R = mult_litt_f<DTYPE, DTYPE>(-tw, ( (In_R[3] - In_R[7]) - (In_I[3] - In_I[7]) ), Q); // CMPX_TWR_TWR_MUL_X_SUB_Y(s18, -0.70710678118654757, X[3], X[7]);
		s18_I = mult_litt_f<DTYPE, DTYPE>(-tw, ( (In_R[3] - In_R[7]) + (In_I[3] - In_I[7]) ), Q);

		t71_R = t66_R + t69_R; // CMPXADD(t71, t66, t69);
		t71_I = t66_I + t69_I;
		t72_R = t66_R - t69_R; // CMPXSUB(t72, t66, t69);
		t72_I = t66_I - t69_I;
		t73_R = t68_R + t70_R; // CMPXADD(t73, t68, t70);
		t73_I = t68_I + t70_I;
		s19_R = t68_I - t70_I; // CMPX_NI_MUL_X_SUB_Y(s19, t68, t70);
		s19_I = t70_R - t68_R;

		Out_R[0] = t71_R + t73_R; // CMPXADD(Y[0], t71, t73);
		Out_I[0] = t71_I + t73_I;
		Out_R[4] = t71_R - t73_R; // CMPXSUB(Y[4], t71, t73);
		Out_I[4] = t71_I - t73_I;
		Out_R[2] = t72_R + s19_R; // CMPXADD(Y[2], t72, s19);
		Out_I[2] = t72_I + s19_I;
		Out_R[6] = t72_R - s19_R; // CMPXSUB(Y[6], t72, s19);
		Out_I[6] = t72_I - s19_I;

		t74_R = t67_R + s17_R; // CMPXADD(t74, t67, s17);
		t74_I = t67_I + s17_I;
		t75_R = t67_R - s17_R; // CMPXSUB(t75, t67, s17);
		t75_I = t67_I - s17_I;
		t76_R = s16_R + s18_R; // CMPXADD(t76, s16, s18);
		t76_I = s16_I + s18_I;
		s20_R = s16_I - s18_I; // CMPX_NI_MUL_X_SUB_Y(s20, s16, s18);
		s20_I = s18_R - s16_R;

		Out_R[1] = t74_R + t76_R; // CMPXADD(Y[1], t74, t76);
		Out_I[1] = t74_I + t76_I;
		Out_R[5] = t74_R - t76_R; // CMPXSUB(Y[5], t74, t76);
		Out_I[5] = t74_I - t76_I;
		Out_R[3] = t75_R + s20_R; // CMPXADD(Y[3], t75, s20);
		Out_I[3] = t75_I + s20_I;
		Out_R[7] = t75_R - s20_R; // CMPXSUB(Y[7], t75, s20);
		Out_I[7] = t75_I - s20_I;
	}
}

/**
 * @brief Data flow utility, converting 2d complex data vector to 1d for un-parallelization purpose
 * 
 * @tparam DTYPE Data Type
 * @tparam SW Streaming width
 * @tparam RADIX FFT Radix
 * @param In_R Input vector (real)
 * @param In_I Input vector (imag)
 * @param Out_R Output vector (real)
 * @param Out_I Output vector (imag)
 */
template<class DTYPE, int SW, int RADIX>
void conv_2to1d(DTYPE In_R[SW/RADIX][RADIX], DTYPE In_I[SW/RADIX][RADIX], DTYPE Out_R[SW], DTYPE Out_I[SW]){
	for(int i = 0; i < SW/RADIX; i++){
		for(int j = 0; j < RADIX; j++){
			Out_R[i*RADIX+j] = In_R[i][j];
			Out_I[i*RADIX+j] = In_I[i][j];
		}
	}
}

/**
 * @brief Data flow utility, converting 1d complex data vector to 2d for parallelization purpose
 * 
 * @tparam DTYPE Data Type
 * @tparam SW Streaming width
 * @tparam RADIX FFT Radix
 * @param In_R Input vector (real)
 * @param In_I Input vector (imag)
 * @param Out_R Output vector (real)
 * @param Out_I Output vector (imag)
 */
template<class DTYPE, int SW, int RADIX>
void conv_1to2d(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE Out_R[SW/RADIX][RADIX], DTYPE Out_I[SW/RADIX][RADIX]){
	for(int i = 0; i < SW/RADIX; i++){
		for(int j = 0; j < RADIX; j++){
			Out_R[i][j] = In_R[i*RADIX+j];
			Out_I[i][j] = In_I[i*RADIX+j];
		}
	}
}

/**
 * @brief Dft stage module, computing all radix 2 dft for SW data
 * 
 * @tparam DTYPE Data Type
 * @tparam SW Streaming Width
 * @param In_R Input vector (real)
 * @param In_I Input vector (imag)
 * @param Out_R Output vector (real)
 * @param Out_I Output vector(imag)
 */
template<class DTYPE, int SW>
void dft_r2(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE Out_R[SW], DTYPE Out_I[SW]){
	#pragma HLS INTERFACE ap_ctrl_none port=return
	#pragma HLS INLINE off
	#pragma HLS PIPELINE

	const unsigned RADIX = 2;
	DTYPE dft_in_R[SW/RADIX][RADIX], dft_in_I[SW/RADIX][RADIX];
	DTYPE dft_out_R[SW/RADIX][RADIX], dft_out_I[SW/RADIX][RADIX];

	conv_1to2d<DTYPE, SW, RADIX>(In_R, In_I, dft_in_R, dft_in_I);
	for(int i = 0; i < SW/RADIX; i++){
		#pragma HLS UNROLL
		radix_core_2(dft_in_R[i], dft_in_I[i], dft_out_R[i], dft_out_I[i]);
	}
	conv_2to1d<DTYPE, SW, RADIX>(dft_out_R, dft_out_I, Out_R, Out_I);
}

/**
 * @brief Dft stage module
 * 
 * @tparam DTYPE Data Type
 * @tparam SW Streaming Width
 * @tparam RADIX Butterfly radix
 * @param In_R Input vector (real)
 * @param In_I Input vector (imag)
 * @param Out_R Output vector (real)
 * @param Out_I Output vector(imag)
 */
template<class DTYPE, int SW, int RADIX>
void dft(DTYPE In_R[SW], DTYPE In_I[SW], DTYPE Out_R[SW], DTYPE Out_I[SW]){
	#pragma HLS INTERFACE ap_ctrl_none port=return
	#pragma HLS INLINE off
	#pragma HLS PIPELINE

	DTYPE dft_in_R[SW/RADIX][RADIX], dft_in_I[SW/RADIX][RADIX];
	DTYPE dft_out_R[SW/RADIX][RADIX], dft_out_I[SW/RADIX][RADIX];

	conv_1to2d<DTYPE, SW, RADIX>(In_R, In_I, dft_in_R, dft_in_I);
	for(int i = 0; i < SW/RADIX; i++){
		#pragma HLS UNROLL
		radix_core<DTYPE, RADIX>(dft_in_R[i], dft_in_I[i], dft_out_R[i], dft_out_I[i]);
	}
	conv_2to1d<DTYPE, SW, RADIX>(dft_out_R, dft_out_I, Out_R, Out_I);
}

/**
 * @brief Spatial permutation module.
 * 
 * @tparam T Data type
 * @tparam SW Streaming width
 * @param spat_in Data to permute
 * @param spat_out Permuted data 
 * @param y_idx Permutation table
 */
template <class T, int SW>
void spatial_permutation(T spat_in[SW], T spat_out[SW], int y_idx[SW]) {
	for (int i=0; i < SW;i++) {
		spat_out[y_idx[i]] = spat_in[i];
	}
}

/**
 * @brief Explicit datapath connection to switch
 * 
 * @tparam T Data type
 * @tparam SW Streaming Width
 * @param con_to_swh_in Input vector
 * @param con_to_swh_out Output vector
 * @param x_idx permutation table
 */
template <class T, int SW>
void connect_to_switch(T con_to_swh_in[SW], T con_to_swh_out[SW/2][2], int x_idx[SW]) {
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_in complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out complete dim=1
	#pragma HLS ARRAY_PARTITION variable=x_idx complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out complete dim=2
	// #pragma HLS ARRAY_PARTITION variable=con_to_swh_out complete

	#pragma HLS INLINE

	for (int i=0; i < SW/2; i++) {
	#pragma HLS UNROLL
		for (int j=0; j < 2; j++) {
			con_to_swh_out[i][j] = con_to_swh_in[x_idx[2*i+j]];
		}
	}
}

/**
 * @brief Explicit datapath connection to switch
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming Width
 * @param con_to_swh_in_R Input real vector
 * @param con_to_swh_in_I Input imag vector
 * @param con_to_swh_out_R Output real vector
 * @param con_to_swh_out_I Output imag vector
 * @param x_idx Permutation table
 */
template <class DTYPE, int SW>
void connect_to_switch(DTYPE con_to_swh_in_R[SW], DTYPE con_to_swh_in_I[SW], DTYPE con_to_swh_out_R[SW/2][2], DTYPE con_to_swh_out_I[SW/2][2], int x_idx[SW]) {
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_in_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_in_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=x_idx complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out_R complete dim=2
	#pragma HLS ARRAY_PARTITION variable=con_to_swh_out_I complete dim=2
	// #pragma HLS ARRAY_PARTITION variable=con_to_swh_out_R complete
	// #pragma HLS ARRAY_PARTITION variable=con_to_swh_out_I complete

	#pragma HLS INLINE

	for (int i=0; i < SW/2; i++) {
	#pragma HLS UNROLL
		for (int j=0; j < 2; j++) {
			con_to_swh_out_R[i][j] = con_to_swh_in_R[x_idx[2*i+j]];
			con_to_swh_out_I[i][j] = con_to_swh_in_I[x_idx[2*i+j]];
		}
	}
}

/**
 * @brief Explicit datapath connection fram switch
 * 
 * @tparam T Data type
 * @tparam SW Streaming Width
 * @param con_from_swh_in Input vector
 * @param con_from_swh_out Output vector
 * @param y_idx Permutation table
 */
template <class T, int SW>
void connect_from_switch(T con_from_swh_in[SW/2][2], T con_from_swh_out[SW], int y_idx[SW]) {
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_out complete dim=1
	#pragma HLS ARRAY_PARTITION variable=y_idx complete dim=1
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in complete dim=2
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_in complete

	#pragma HLS INLINE

	for (int i=0; i < SW/2; i++) {
	#pragma HLS UNROLL
		for (int j=0; j < 2; j++) {
			con_from_swh_out[y_idx[2*i+j]] = con_from_swh_in[i][j];
		}
	}
}

/**
 * @brief Explicit datapath connection from switch
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming Width
 * @param con_from_swh_in_R Input real Vector
 * @param con_from_swh_in_I Input imag Vector
 * @param con_from_swh_out_R Output real Vector
 * @param con_from_swh_out_I Output imag Vector
 * @param y_idx Permutation table
 */
template <class DTYPE, int SW>
void connect_from_switch(DTYPE con_from_swh_in_R[SW/2][2], DTYPE con_from_swh_in_I[SW/2][2], DTYPE con_from_swh_out_R[SW],  DTYPE con_from_swh_out_I[SW], int y_idx[SW]) {
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in_R complete dim=1
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_out_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_out_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=y_idx complete dim=1
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in_R complete dim=2
	// #pragma HLS ARRAY_PARTITION variable=con_from_swh_in_I complete dim=2
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_in_R complete
	#pragma HLS ARRAY_PARTITION variable=con_from_swh_in_I complete

	#pragma HLS INLINE

	for (int i=0; i < SW/2; i++) {
	#pragma HLS UNROLL
		for (int j=0; j < 2; j++) {
			con_from_swh_out_R[y_idx[2*i+j]] = con_from_swh_in_R[i][j];
			con_from_swh_out_I[y_idx[2*i+j]] = con_from_swh_in_I[i][j];
		}
	}
}

/**
 * @brief Switch network on two data
 * 
 * @tparam T Data type
 * @param in Input vector
 * @param out Output vector
 * @param switch_on Switch enable
 */
template <class T>
void switch_2_by_2(T in[2], T out[2], bool switch_on) {
	out[0] = in[(switch_on==true)?(1):(0)];
	out[1] = in[(switch_on==true)?(0):(1)];
}

/**
 * @brief Switch network on two data
 * 
 * @tparam DTYPE Data type
 * @param in_R Input real vector
 * @param in_I Input imag vector
 * @param out_R Output real vector
 * @param out_I Output imag vector
 * @param switch_on Switch enable
 */
template <class DTYPE>
void switch_2_by_2(DTYPE in_R[2], DTYPE in_I[2], DTYPE out_R[2], DTYPE out_I[2], bool switch_on) {
	out_R[0] = in_R[(switch_on==true)?(1):(0)];
	out_I[0] = in_I[(switch_on==true)?(1):(0)];
	out_R[1] = in_R[(switch_on==true)?(0):(1)];
	out_I[1] = in_I[(switch_on==true)?(0):(1)];
}

/**
 * @brief Switch array on SW data simultaneously
 * 
 * @tparam T Data type
 * @tparam SW Streaming Width
 * @param num amount of switches
 * @param swh_array_in Input vector
 * @param swh_array_out Output vector
 * @param switch_on Switch enable
 */
template <class T, int SW>
void switch_array(int num, T swh_array_in[SW/2][2], T swh_array_out[SW/2][2], bool switch_on) {
	for (int i=0; i<num; i++) {
	#pragma HLS UNROLL
		switch_2_by_2<T>(swh_array_in[i], swh_array_out[i], switch_on);
	}
}

/**
 * @brief Switch array on SW data simultaneously
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming Width
 * @param num amount of switches
 * @param swh_array_in_R Input real vector
 * @param swh_array_in_I Input imag vector
 * @param swh_array_out_R Output real vector
 * @param swh_array_out_I Output imag vector
 * @param switch_on Switch enable
 */
template <class DTYPE, int SW>
void switch_array(int num, DTYPE swh_array_in_R[SW/2][2], DTYPE swh_array_in_I[SW/2][2], DTYPE swh_array_out_R[SW/2][2], DTYPE swh_array_out_I[SW/2][2], bool switch_on) {
	for (int i=0; i<num; i++) {
	#pragma HLS UNROLL
		switch_2_by_2<DTYPE>(swh_array_in_R[i], swh_array_in_I[i], swh_array_out_R[i],swh_array_out_I[i], switch_on);
	}
}

/**
 * @brief Stage of switch
 * 
 * @tparam T Data type
 * @tparam SW Streaming width
 * @param PRE_IN Input vector
 * @param NEXT_PRE_IN Output vector
 * @param control Switch enable
 * @param idx_in_or_out Permutation table
 */
template <class T, int SW>
void onestage(T PRE_IN[SW], T NEXT_PRE_IN[SW], bool control, int idx_in_or_out[SW]) {
	#pragma HLS INLINE
    T IN[SW/2][2], OUT[SW/2][2];

    connect_to_switch<T, SW>(PRE_IN, IN, idx_in_or_out);
    switch_array<T, SW>(SW/2, IN, OUT, control);
    connect_from_switch<T, SW>(OUT, NEXT_PRE_IN, idx_in_or_out);
}

/**
 * @brief Stage of switch
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming width
 * @param PRE_IN_R Input real vector
 * @param PRE_IN_I Input imag vector
 * @param NEXT_PRE_IN_R Output real vector
 * @param NEXT_PRE_IN_I Output imag vector
 * @param control Switch enable
 * @param idx_in_or_out Permutation table
 */
template <class DTYPE, int SW>
void onestage(DTYPE PRE_IN_R[SW], DTYPE PRE_IN_I[SW], DTYPE NEXT_PRE_IN_R[SW], DTYPE NEXT_PRE_IN_I[SW], bool control, int idx_in_or_out[SW]) {
	#pragma HLS INLINE
    DTYPE IN_R[SW/2][2], IN_I[SW/2][2], OUT_R[SW/2][2], OUT_I[SW/2][2];

    connect_to_switch<DTYPE, SW>(PRE_IN_R, PRE_IN_I, IN_R, IN_I, idx_in_or_out);
    switch_array<DTYPE, SW>(SW/2, IN_R, IN_I, OUT_R, OUT_I, control);
    connect_from_switch<DTYPE, SW>(OUT_R, OUT_I, NEXT_PRE_IN_R, NEXT_PRE_IN_I, idx_in_or_out);
}

/**
 * @brief Vector copy helper
 * 
 * @tparam T Data type
 * @tparam SW Streaming Width
 * @param from Input vector
 * @param to Output vector
 */
template <class T, int SW>
void copy_vector(T from[SW], T to[SW]) {
	for (int i=0; i < SW; i++) {
		to[i] = from[i];
	}
}

/**
 * @brief Vector copy helper
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming width
 * @param from_R Input real vector
 * @param from_I Input imag vector
 * @param to_R Output real vector
 * @param to_I Output imag vector
 */
template <class DTYPE, int SW>
void copy_vector(DTYPE from_R[SW], DTYPE from_I[SW], DTYPE to_R[SW], DTYPE to_I[SW]) {
	for (int i=0; i < SW; i++) {
		to_R[i] = from_R[i];
		to_I[i] = from_I[i];
	}
}

/**
 * @brief Switch Network Input function
 * 
 * @tparam T Data type
 * @tparam NUM_STAGE Stage ID
 * @tparam SW Streaming Width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 * @param X Input vector
 * @param Y Output vector
 * @param bits 
 * @param init_perm_idx 
 * @param idx_in_or_out 
 * @param control_bit 
 */
template <class T, int NUM_STAGE, int SW, int LOG2N, int LOG2SW>
void switch_network_write(T X[SW], T Y[SW], ap_uint<LOG2N-LOG2SW> bits, int init_perm_idx[SW], int idx_in_or_out[NUM_STAGE][SW], int control_bit[NUM_STAGE]) {
	#pragma HLS ARRAY_PARTITION variable=X complete dim=1
	#pragma HLS ARRAY_PARTITION variable=Y complete dim=1
	#pragma HLS INTERFACE ap_ctrl_none port=return
	#pragma HLS INLINE off

	T PRE_IN[NUM_STAGE+1][SW];
	#pragma HLS ARRAY_PARTITION variable=PRE_IN complete dim=1

	// initial permutation
	spatial_permutation<T, SW>(X, PRE_IN[0], init_perm_idx);

	for (int i=0; i < NUM_STAGE; i++) {
	#pragma HLS UNROLL
	    onestage<T, SW>(PRE_IN[i], PRE_IN[i+1], bits.get_bit(control_bit[i]-LOG2SW), idx_in_or_out[i]);
	}

	// output port
	copy_vector<T, SW>(PRE_IN[NUM_STAGE], Y);
}

/**
 * @brief Switch Network output function
 * 
 * @tparam DTYPE Data type
 * @tparam NUM_STAGE Stage ID
 * @tparam SW Streaming Width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 * @param X_R Input real Vector
 * @param X_I Input imag Vector
 * @param Y_R Output real Vector
 * @param Y_I Output imag Vector
 * @param bits 
 * @param idx_in_or_out 
 * @param control_bit 
 */
template <typename DTYPE, int NUM_STAGE, int SW, int LOG2N, int LOG2SW>
void switch_network_read(DTYPE X_R[SW], DTYPE X_I[SW], DTYPE Y_R[SW], DTYPE Y_I[SW], ap_uint<LOG2N-LOG2SW> bits, int idx_in_or_out[NUM_STAGE][SW], int control_bit[NUM_STAGE]) {
	#pragma HLS ARRAY_PARTITION variable=X_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=X_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=Y_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=Y_I complete dim=1
	#pragma HLS INTERFACE ap_ctrl_none port=return
	#pragma HLS INLINE off

	DTYPE PRE_IN_R[NUM_STAGE+1][SW], PRE_IN_I[NUM_STAGE+1][SW];
	#pragma HLS ARRAY_PARTITION variable=PRE_IN_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=PRE_IN_I complete dim=1

	// no initial permutation for the read switch network.
	copy_vector<DTYPE, SW>(X_R, X_I, PRE_IN_R[0], PRE_IN_I[0]);

	for (int i=0; i<NUM_STAGE; i++) {
		#pragma HLS UNROLL
	    onestage<DTYPE, SW>(PRE_IN_R[i], PRE_IN_I[i], PRE_IN_R[i+1], PRE_IN_I[i+1], bits.get_bit(control_bit[i]-LOG2SW), idx_in_or_out[i]);
	}

	// output port
	copy_vector<DTYPE, SW>(PRE_IN_R[NUM_STAGE], PRE_IN_I[NUM_STAGE], Y_R, Y_I);
}

/**
 * @brief Buffer write address generator
 * 
 * @tparam SW Streaming Width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 * @param in_count 
 * @param flip 
 * @param bit_seq 
 * @param widx 
 */
template<int SW, int LOG2N, int LOG2SW>
void buf_write_addr_generation(ap_uint<LOG2N-LOG2SW> in_count, bool flip, int bit_seq[LOG2N-LOG2SW], ap_uint<LOG2N-LOG2SW+1> widx[2]) {
#pragma HLS INLINE
	int i, j;
	for (i=0; i<SW; i++) {
	#pragma HLS UNROLL
		ap_uint<LOG2N> bits = SW*in_count+i;

		for (j=0; j<LOG2N-LOG2SW; j++) {
	#pragma HLS UNROLL
			bool bit = bits.get_bit( bit_seq[LOG2N-LOG2SW-1-j] );
			widx[i].set_bit(j, bit);
		}
		widx[i].set_bit(LOG2N-LOG2SW, flip);
	}
}

/**
 * @brief Buffer read address generator
 * 
 * @tparam SW Streaming Width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 * @param out_count 
 * @param flip 
 * @param ridx 
 */
template<int SW, int LOG2N, int LOG2SW>
void buf_read_addr_generation(ap_uint<LOG2N-LOG2SW> out_count, bool flip, ap_uint<LOG2N-LOG2SW+1> ridx[2]) {
	int i;
	for (i=0; i<SW; i++) {
		#pragma HLS UNROLL
		ridx[i].range(LOG2N-LOG2SW-1, 0) = out_count.range(LOG2N-LOG2SW-1, 0);
		ridx[i].set_bit(LOG2N-LOG2SW, flip);
	}
}

/**
 * @brief Address/Data combiner helper
 * 
 * @tparam DTYPE Data type
 * @tparam SW Streaming Width
 * @tparam LOG2N log2(N)
 * @tparam LOG2SW log2(SW)
 * @param combination Output
 * @param addr adress
 * @param data_R Input real vector
 * @param data_I Input imag vector
 */
template<class DTYPE, int SW, int LOG2N, int LOG2SW>
void combine_addr_data(content_addr<DTYPE, LOG2N, LOG2SW> combination[SW], ap_uint<LOG2N-LOG2SW+1> addr[SW], DTYPE data_R[SW], DTYPE data_I[SW]) {
	int i;
	for (i = 0; i < SW; i++) {
		#pragma HLS UNROLL
		combination[i].data_R = data_R[i];
		combination[i].data_I = data_I[i];
		combination[i].addr = addr[i];
	}
}

#endif //SPIRAL_UTILS_HPP_
