/**
 * @file lib_fft_sbg_radix.hpp
 * @author Hugues Almorin (hugues.almorin@arelis.com)
 * @brief This file contains software wrapping for using the FFT model
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef LIB_FFT_SBG_RADIX_HPP_
#define LIB_FFT_SBG_RADIX_HPP_

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include "../../common/Parameters.hpp"
#include "../../common/DataVector.hpp"
#include "../module/fft_sbg_radix.hpp"
#include "ap_int.h"


void lib_fft_sbg_radix_load_fft_roms(std::string cfile, const int32_t SIZE);
void lib_sbg_radix_fft(DataVector& data, Parameters& p);

#endif // LIB_FFT_SBG_RADIX_HPP_