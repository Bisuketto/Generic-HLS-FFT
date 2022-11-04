/**
 * @file fft_sbg_radix_hlstop.hpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains the HLS top module for the FFT
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef FFT_SBG_RADIX_HLSTOP_HPP_
#define FFT_SBG_RADIX_HLSTOP_HPP_

#include "ap_int.h"
#include "fft_sbg_radix.hpp"

#ifndef NFFT
#define NFFT 		4096
#endif

#ifndef TW_WIDTH
#define TW_WIDTH    20
#endif

#ifndef STREAMING_WIDTH
#define STREAMING_WIDTH 8
#endif

#ifndef GROUP_SIZE
#define GROUP_SIZE 2
#endif

#ifndef CORE_RADIX
#define CORE_RADIX 8
#endif

#ifndef COS_LUT
#define COS_LUT_CONCAT(N, I, F) COS_LUT_N##N##_Q##I##_##F
#define COS_LUT_HELPER(N, I, F) COS_LUT_CONCAT(N, I, F)
#define COS_LUT COS_LUT_N4096_Q2_18 //COS_LUT_HELPER(NFFT, 2, TW_WIDTH-2)
#endif

#ifndef SIN_LUT
#define SIN_LUT_CONCAT(N, I, F) SIN_LUT_N##N##_Q##I##_##F
#define SIN_LUT_HELPER(N, I, F) SIN_LUT_CONCAT(N, I, F)
#define SIN_LUT SIN_LUT_N4096_Q2_18 //SIN_LUT_HELPER(NFFT, 2, TW_WIDTH-2)
#endif

#ifndef DIGREV_CONFIG
#warning "Using default DIGREV_CONFIG"
#define DIGREV_CONFIG_CONCAT(N, SW, R) DIG_REV_PERM_CONFIG_N##N##_SW##SW##_R##R
#define DIGREV_CONFIG_HELPER(N, SW, R) DIGREV_CONFIG_CONCAT(N, SW, R)
#define DIGREV_CONFIG DIGREV_CONFIG_HELPER(NFFT, STREAMING_WIDTH, CORE_RADIX)
#endif

#ifndef DIGREV_NUMSTAGE
#warning "Using default DIGREV_NUMSTAGE"
#define DIGREV_NUMSTAGE_CONCAT(N, SW, R) DIGIT_REV_NUM_STAGE_N##N##_SW##SW##_R##R
#define DIGREV_NUMSTAGE_HELPER(N, SW, R) DIGREV_NUMSTAGE_CONCAT(N, SW, R)
#define DIGREV_NUMSTAGE DIGREV_NUMSTAGE_HELPER(NFFT, STREAMING_WIDTH, CORE_RADIX)
#endif

#ifndef STRIDE_CONFIG
#warning "Using default STRIDE_CONFIG"
#define STRIDE_CONFIG_CONCAT(N, SW, R) STRIDE_PERM_CONFIG_N##N##_SW##SW##_R##R
#define STRIDE_CONFIG_HELPER(N, SW, R) STRIDE_CONFIG_CONCAT(N, SW, R)
#define STRIDE_CONFIG STRIDE_CONFIG_HELPER(NFFT, STREAMING_WIDTH, CORE_RADIX)
#endif

#ifndef STRIDE_NUMSTAGE
#warning "Using default STRIDE_NUMSTAGE"
#define STRIDE_NUMSTAGE_CONCAT(N, SW, R) STRIDE_PERM_SWITCH_NUM_STAGE_N##N##_SW##SW##_R##R
#define STRIDE_NUMSTAGE_HELPER(N, SW, R) STRIDE_NUMSTAGE_CONCAT(N, SW, R)
#define STRIDE_NUMSTAGE STRIDE_NUMSTAGE_HELPER(NFFT, STREAMING_WIDTH, CORE_RADIX)
#endif

#define D_WIDTH		22
#define I_WIDTH		TW_WIDTH+D_WIDTH

typedef ap_int<I_WIDTH> SAMPLE_TYPE;
typedef ap_int<TW_WIDTH> TW_TYPE;
typedef ap_int<I_WIDTH> INT_TYPE;

void fft(SAMPLE_TYPE din_I[NFFT], SAMPLE_TYPE din_Q[NFFT], SAMPLE_TYPE dout_I[NFFT], SAMPLE_TYPE dout_Q[NFFT]);

#endif // FFT_SBG_RADIX_HLSTOP_HPP_
