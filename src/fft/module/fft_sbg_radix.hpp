/**
 * @file fft_sbg_radix.hpp
 * @author Hugues Almorin (hugues.almorin@arelis.com)
 * @brief This file contains the HLS-ready FFT description allowing Butterfly and Stage Group parallelization, multi radix. Based on Spiral works
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */


#ifndef FFT_SBG_RADIX_HPP_
#define FFT_SBG_RADIX_HPP_

#include "../../common/const_logs.hpp"
#include "../../../roms/roms.h"
#include "fft_sbg_radix_utils.hpp"
#include "ap_int.h"


template<class DTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GSTART, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_st_core(	DTYPE buf_in_R[SW][SIZE/SW], DTYPE buf_in_I[SW][SIZE/SW], DTYPE buf_out_R[SW][SIZE/SW], DTYPE buf_out_I[SW][SIZE/SW], TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config, int STAGE, int Q){
	#pragma HLS INLINE
	#pragma HLS ARRAY_PARTITION variable=buf_in_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_in_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_out_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_out_I complete dim=1

	const unsigned i_stage = STAGE;

	DTYPE din_R[SW], din_I[SW];
	DTYPE dout_R[SW], dout_I[SW];
	DTYPE twiddled_data_R[SW], twiddled_data_I[SW];

	const bool rd_flip = false;
	const bool wt_flip = false;
	const bool first_stage = (i_stage+GSTART == 0)?(true):(false);

	butterfly_loop :
	for(int j = 0; j < SIZE/SW; j++){
		#pragma HLS PIPELINE
		fft_sbg_radix_st_dr_buf_read<DTYPE, SIZE, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, rd_flip, first_stage, din_R, din_I, dig_rev_config, stride_config);


		fft_sbg_radix_twiddling<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART>(din_R, din_I, twiddled_data_R, twiddled_data_I, Tw_R, Tw_I, i_stage+1, Q);
		dft<DTYPE, SW, RADIX>(twiddled_data_R, twiddled_data_I, dout_R, dout_I);

		fft_sbg_radix_stride_buf_write<DTYPE, SIZE, SW, GSTART, LOG2N, LOG2SW, STRIDE_PERM_SWITCH_NUM_STAGE>(dout_R, dout_I, buf_out_R, buf_out_I, wt_flip, stride_config);
	}	
}

template<class DTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GROUP, int GSTART, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_grp(DTYPE buf_in_R[SW][SIZE/SW], DTYPE buf_in_I[SW][SIZE/SW], DTYPE buf_out_R[SW][SIZE/SW], DTYPE buf_out_I[SW][SIZE/SW],
			TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config, int Q){
	
	static DTYPE buf_inter_R_0[SW][SIZE/SW], buf_inter_I_0[SW][SIZE/SW];
	static DTYPE buf_inter_R_1[SW][SIZE/SW], buf_inter_I_1[SW][SIZE/SW];
	static DTYPE buf_inter_R_2[SW][SIZE/SW], buf_inter_I_2[SW][SIZE/SW];
	static DTYPE buf_inter_R_3[SW][SIZE/SW], buf_inter_I_3[SW][SIZE/SW];
	static DTYPE buf_inter_R_4[SW][SIZE/SW], buf_inter_I_4[SW][SIZE/SW];
	static DTYPE buf_inter_R_5[SW][SIZE/SW], buf_inter_I_5[SW][SIZE/SW];
	static DTYPE buf_inter_R_6[SW][SIZE/SW], buf_inter_I_6[SW][SIZE/SW];
	static DTYPE buf_inter_R_7[SW][SIZE/SW], buf_inter_I_7[SW][SIZE/SW];
	static DTYPE buf_inter_R_8[SW][SIZE/SW], buf_inter_I_8[SW][SIZE/SW];
	static DTYPE buf_inter_R_9[SW][SIZE/SW], buf_inter_I_9[SW][SIZE/SW];
	static DTYPE buf_inter_R_10[SW][SIZE/SW], buf_inter_I_10[SW][SIZE/SW];

	switch(GROUP){
		case 1:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			break;
		case 2:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			break;
		case 3:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			break;
		case 4:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			break;
		case 5:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			break;
		case 6:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			break;
		case 7:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			break;
		case 8:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_inter_R_6, buf_inter_I_6, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_6, buf_inter_I_6, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+7, Q);
			break;
		case 9:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_inter_R_6, buf_inter_I_6, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_6, buf_inter_I_6, buf_inter_R_7, buf_inter_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+7, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_7, buf_inter_I_7, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+8, Q);
			break;
		case 10:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_inter_R_6, buf_inter_I_6, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_6, buf_inter_I_6, buf_inter_R_7, buf_inter_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+7, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_7, buf_inter_I_7, buf_inter_R_8, buf_inter_I_8, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+8, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_8, buf_inter_I_8, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+9, Q);
			break;
		case 11:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_inter_R_6, buf_inter_I_6, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_6, buf_inter_I_6, buf_inter_R_7, buf_inter_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+7, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_7, buf_inter_I_7, buf_inter_R_8, buf_inter_I_8, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+8, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_8, buf_inter_I_8, buf_inter_R_9, buf_inter_I_9, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+9, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_9, buf_inter_I_9, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+10, Q);
			break;
		case 12:
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R_0, buf_inter_I_0, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+0, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_0, buf_inter_I_0, buf_inter_R_1, buf_inter_I_1, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+1, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_1, buf_inter_I_1, buf_inter_R_2, buf_inter_I_2, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+2, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_2, buf_inter_I_2, buf_inter_R_3, buf_inter_I_3, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+3, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_3, buf_inter_I_3, buf_inter_R_4, buf_inter_I_4, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+4, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_4, buf_inter_I_4, buf_inter_R_5, buf_inter_I_5, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+5, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_5, buf_inter_I_5, buf_inter_R_6, buf_inter_I_6, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+6, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_6, buf_inter_I_6, buf_inter_R_7, buf_inter_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+7, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_7, buf_inter_I_7, buf_inter_R_8, buf_inter_I_8, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+8, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_8, buf_inter_I_8, buf_inter_R_9, buf_inter_I_9, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+9, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_9, buf_inter_I_9, buf_inter_R_10, buf_inter_I_10, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+10, Q);
			sbg_radix_st_core<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GSTART, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R_10, buf_inter_I_10, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, GSTART+11, Q);
			break;
	}
}

template<class DTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GROUP, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_stages(DTYPE buf_in_R[SW][SIZE/SW], DTYPE buf_in_I[SW][SIZE/SW], DTYPE buf_out_R[SW][SIZE/SW], DTYPE buf_out_I[SW][SIZE/SW],
			TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config, int Q){

	#pragma HLS INLINE
	
	constexpr const unsigned M = ceillogR(SIZE, RADIX);
	constexpr const unsigned divi = M % GROUP;
	constexpr const unsigned intstages = ((M/GROUP >= 1) && (divi == 0)) ? M/GROUP : M/GROUP+1;
	constexpr const unsigned n_inter = (intstages > 1) ? intstages - 1 : 1; // avoiding Inter_X dim_1 = 0
	static_assert(M >= GROUP, "GROUP SHOULD BE <= N_STAGES");

	static DTYPE buf_inter_R[n_inter][SW][SIZE/SW], buf_inter_I[n_inter][SW][SIZE/SW];
	#pragma HLS ARRAY_PARTITION variable=buf_inter_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_inter_I complete dim=1

	if(divi == 0){
		switch(intstages){
			case 1:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 2:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 3:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 4:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 5:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_inter_R[3], buf_inter_I[3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[3], buf_inter_I[3], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 6:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_inter_R[3], buf_inter_I[3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[3], buf_inter_I[3], buf_inter_R[4], buf_inter_I[4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[4], buf_inter_I[4], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 7:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 8:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 9:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 10:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 11:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_inter_R[ 9], buf_inter_I[ 9], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP,10*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 9], buf_inter_I[ 9], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 12:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_inter_R[ 9], buf_inter_I[ 9], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP,10*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 9], buf_inter_I[ 9], buf_inter_R[10], buf_inter_I[10], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP,11*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[10], buf_inter_I[10], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
		}
	}
	else{
		switch(intstages){
			case 1:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 2:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 3:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 4:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 5:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_inter_R[3], buf_inter_I[3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[3], buf_inter_I[3], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 6:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[0], buf_inter_I[0], buf_inter_R[1], buf_inter_I[1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[1], buf_inter_I[1], buf_inter_R[2], buf_inter_I[2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[2], buf_inter_I[2], buf_inter_R[3], buf_inter_I[3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[3], buf_inter_I[3], buf_inter_R[4], buf_inter_I[4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[4], buf_inter_I[4], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 7:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 8:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 9:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 10:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 11:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_inter_R[ 9], buf_inter_I[ 9], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi,10*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 9], buf_inter_I[ 9], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
			case 12:
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 0*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_in_R, buf_in_I, buf_inter_R[0], buf_inter_I[0], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 1*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 0], buf_inter_I[ 0], buf_inter_R[ 1], buf_inter_I[ 1], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 2*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 1], buf_inter_I[ 1], buf_inter_R[ 2], buf_inter_I[ 2], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 3*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 2], buf_inter_I[ 2], buf_inter_R[ 3], buf_inter_I[ 3], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 4*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 3], buf_inter_I[ 3], buf_inter_R[ 4], buf_inter_I[ 4], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 5*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 4], buf_inter_I[ 4], buf_inter_R[ 5], buf_inter_I[ 5], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 6*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 5], buf_inter_I[ 5], buf_inter_R[ 6], buf_inter_I[ 6], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 7*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 6], buf_inter_I[ 6], buf_inter_R[ 7], buf_inter_I[ 7], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 8*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 7], buf_inter_I[ 7], buf_inter_R[ 8], buf_inter_I[ 8], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, 9*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 8], buf_inter_I[ 8], buf_inter_R[ 9], buf_inter_I[ 9], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP,10*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[ 9], buf_inter_I[ 9], buf_inter_R[10], buf_inter_I[10], Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				sbg_radix_grp<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW,  divi,11*GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_inter_R[10], buf_inter_I[10], buf_out_R, buf_out_I, Tw_R, Tw_I, dig_rev_config, stride_config, Q);
				break;
		}
	}
}

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE>
void sbg_radix_buff_fill(DTYPE In_R[SIZE], DTYPE In_I[SIZE], DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config){
	DTYPE packeted_in_R[SW], packeted_in_I[SW];
	#pragma HLS ARRAY_PARTITION variable=packeted_in_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=packeted_in_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_I complete dim=1
	buf_fill_loop : 
	for(int i = 0; i < SIZE/SW; i++){
		#pragma HLS PIPELINE
		read_in_loop :
		for(int j = 0; j < SW; j++){
			packeted_in_R[j] = In_R[SW*i+j];
			packeted_in_I[j] = In_I[SW*i+j];
		}
		fft_sbg_radix_digrev_buf_write<DTYPE, SIZE, SW, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE>(packeted_in_R, packeted_in_I, buf_R, buf_I, false, dig_rev_config);
	}
}

template<class DTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_buff_dump(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], DTYPE Out_R[SIZE], DTYPE Out_I[SIZE], perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config){
	DTYPE packeted_out_R[SW], packeted_out_I[SW];
	#pragma HLS ARRAY_PARTITION variable=packeted_out_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=packeted_out_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_I complete dim=1
	constexpr const bool final_rd_flip = false;
	buf_dump_loop :
	for(int i = 0; i < SIZE/SW; i++){
		#pragma HLS PIPELINE
		fft_sbg_radix_final_stride_buf_read<DTYPE, SIZE, SW, LOG2N, LOG2SW, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R, buf_I, final_rd_flip, packeted_out_R, packeted_out_I, stride_config);

		write_out_loop :
		for(int j = 0; j < SW; j++){
			Out_R[SW*i+j] = packeted_out_R[j];
			Out_I[SW*i+j] = packeted_out_I[j];
		}
	}
}

template<class DTYPE, class STYPE, class VTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE>
void sbg_radix_buff_fill_axis(STYPE &In, DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config){
	DTYPE packeted_in_R[SW], packeted_in_I[SW];
	VTYPE in_vector;
	#pragma HLS ARRAY_PARTITION variable=packeted_in_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=packeted_in_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_I complete dim=1
	buf_fill_loop : 
	for(int i = 0; i < SIZE/SW; i++){
		#pragma HLS PIPELINE
		read_in_loop :
		for(int j = 0; j < SW; j++){
			In.read(in_vector);

			// packeted_in_R[j] = In_R[SW*i+j];
			packeted_in_R[j] = (in_vector.data).range(31, 0);
			// packeted_in_I[j] = In_I[SW*i+j];
			packeted_in_I[j] = (in_vector.data).range(63, 32);
		}
		fft_sbg_radix_digrev_buf_write<DTYPE, SIZE, SW, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE>(packeted_in_R, packeted_in_I, buf_R, buf_I, false, dig_rev_config);
	}
}

template<class DTYPE, class STYPE, class VTYPE, int SIZE, int SW, int LOG2N, int LOG2SW, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_buff_dump_axis(DTYPE buf_R[SW][SIZE/SW], DTYPE buf_I[SW][SIZE/SW], STYPE &Out, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config){
	DTYPE packeted_out_R[SW], packeted_out_I[SW];
	VTYPE out_vector;
	#pragma HLS ARRAY_PARTITION variable=packeted_out_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=packeted_out_I complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_R complete dim=1
	#pragma HLS ARRAY_PARTITION variable=buf_I complete dim=1
	constexpr const bool final_rd_flip = false;
	buf_dump_loop :
	for(int i = 0; i < SIZE/SW; i++){
		#pragma HLS PIPELINE
		fft_sbg_radix_final_stride_buf_read<DTYPE, SIZE, SW, LOG2N, LOG2SW, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R, buf_I, final_rd_flip, packeted_out_R, packeted_out_I, stride_config);

		write_out_loop :
		for(int j = 0; j < SW; j++){
			// Out_R[SW*i+j] = packeted_out_R[j];
			// Out_I[SW*i+j] = packeted_out_I[j];
			out_vector.data = (packeted_out_I[j],packeted_out_R[j]);
			out_vector.keep = -1;
			out_vector.strb = -1;
			out_vector.last = ((SW*i+j) >= (SIZE-1)) ? 1 : 0;
			Out.write(out_vector);
		}
	}
}

/**
 * @brief FFT, using permutations and explicit datapath
 * 
 * @tparam DTYPE Data type
 * @tparam TTYPE Twiddle factor type
 * @tparam SIZE FFT Size
 * @tparam RADIX Butterfly radix
 * @tparam SW Streaming Width
 * @param In_R Input vector (real)
 * @param In_I Input vector (imag)
 * @param Out_R Output vector (real)
 * @param Out_I Output vector (imag)
 * @param Tw_R Twiddle factor (real)
 * @param Tw_I Twiddle factor (imag)
 */
template<class DTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GROUP, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_fft(DTYPE In_R[SIZE], DTYPE In_I[SIZE], DTYPE Out_R[SIZE], DTYPE Out_I[SIZE], TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config, int Q){
	#pragma HLS INLINE

	#pragma HLS ARRAY_PARTITION variable=In_R cyclic factor=SW dim=1
	#pragma HLS ARRAY_PARTITION variable=In_I cyclic factor=SW dim=1
	#pragma HLS ARRAY_PARTITION variable=Out_R cyclic factor=SW dim=1
	#pragma HLS ARRAY_PARTITION variable=Out_I cyclic factor=SW dim=1

	static DTYPE buf_R_0[SW][SIZE/SW], buf_I_0[SW][SIZE/SW];
	static DTYPE buf_R_7[SW][SIZE/SW], buf_I_7[SW][SIZE/SW];

	sbg_radix_buff_fill<DTYPE, SIZE, SW, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE>(In_R, In_I, buf_R_0, buf_I_0, dig_rev_config);

	sbg_radix_stages<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R_0, buf_I_0, buf_R_7, buf_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, Q);

	sbg_radix_buff_dump<DTYPE, SIZE, SW, LOG2N, LOG2SW, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R_7, buf_I_7, Out_R, Out_I, stride_config);

}

/**
 * @brief FFT, using permutations and explicit datapath with axis interfaces
 * 
 * @tparam DTYPE data type
 * @tparam STYPE stream type
 * @tparam VTYPE vector type
 * @tparam TTYPE Twiddle factor type
 * @tparam SIZE FFT Size
 * @tparam RADIX Butterfly radix
 * @tparam SW Streaming Width
 * @param In Input stream (real;imag)
 * @param Out Output stream (real; imag)
 * @param Tw_R Twiddle factor (real)
 * @param Tw_I Twiddle factor (imag)
 */
template<class DTYPE, class STYPE, class VTYPE, class TTYPE, class ITYPE, int SIZE, int RADIX, int SW, int GROUP, int LOG2N, int LOG2SW, int DIGIT_REV_NUM_STAGE, int STRIDE_PERM_SWITCH_NUM_STAGE>
void sbg_radix_fft_axis(STYPE &din, STYPE &dout, TTYPE Tw_R[SIZE], TTYPE Tw_I[SIZE], perm_config<DIGIT_REV_NUM_STAGE, SW, LOG2N, LOG2SW> dig_rev_config, perm_config<STRIDE_PERM_SWITCH_NUM_STAGE, SW, LOG2N, LOG2SW> stride_config, int Q){
	#pragma HLS INLINE

	static DTYPE buf_R_0[SW][SIZE/SW], buf_I_0[SW][SIZE/SW];
	static DTYPE buf_R_7[SW][SIZE/SW], buf_I_7[SW][SIZE/SW];

	sbg_radix_buff_fill_axis<DTYPE, STYPE, VTYPE, SIZE, SW, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE>(din, buf_R_0, buf_I_0, dig_rev_config);

	sbg_radix_stages<DTYPE, TTYPE, ITYPE, SIZE, RADIX, SW, GROUP, LOG2N, LOG2SW, DIGIT_REV_NUM_STAGE, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R_0, buf_I_0, buf_R_7, buf_I_7, Tw_R, Tw_I, dig_rev_config, stride_config, Q);

	sbg_radix_buff_dump_axis<DTYPE, STYPE, VTYPE, SIZE, SW, LOG2N, LOG2SW, STRIDE_PERM_SWITCH_NUM_STAGE>(buf_R_7, buf_I_7, dout, stride_config);

}

#endif // FFT_SBG_RADIX_HPP_
