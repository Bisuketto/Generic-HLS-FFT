/**
 * @file fft_sbg_radix_hlstop.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains the HLS top module for the FFT
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include "fft_sbg_radix_hlstop.hpp"

void fft(SAMPLE_TYPE din_I[NFFT], SAMPLE_TYPE din_Q[NFFT], SAMPLE_TYPE dout_I[NFFT], SAMPLE_TYPE dout_Q[NFFT]){
	#pragma HLS DATAFLOW


	static TW_TYPE t_cos_fixed[NFFT] = COS_LUT;

	static TW_TYPE t_sin_fixed[NFFT] = SIN_LUT;

	static perm_config<DIGREV_NUMSTAGE, STREAMING_WIDTH, ceillog2(NFFT), ceillog2(STREAMING_WIDTH)> digrev_config = {DIGREV_CONFIG};

	static perm_config<STRIDE_NUMSTAGE, STREAMING_WIDTH, ceillog2(NFFT), ceillog2(STREAMING_WIDTH)> stride_config = {STRIDE_CONFIG};

	sbg_radix_fft<SAMPLE_TYPE, TW_TYPE, INT_TYPE,  NFFT, CORE_RADIX, STREAMING_WIDTH, GROUP_SIZE, ceillog2(NFFT), ceillog2(STREAMING_WIDTH), DIGREV_NUMSTAGE , STRIDE_NUMSTAGE >(din_I, din_Q, dout_I, dout_Q, t_cos_fixed, t_sin_fixed, digrev_config, stride_config, TW_WIDTH-3);
}
