# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(fft_fftw LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(fft_fftw
    sw_wrap/lib_fft_float_gold_fftw.cpp
)
target_link_libraries(fft_fftw fftw3 fft_tools)
target_include_directories( fft_fftw PRIVATE
        ${AP_TYPES_LIB_DIR}
)