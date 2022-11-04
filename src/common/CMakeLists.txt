# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(fft_tools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(fft_tools 
    DataVector.cpp
    mult_add.cpp
    Parameters.cpp
    fft_fftw3.cpp
)
target_include_directories( fft_tools PRIVATE 
    ${AP_TYPES_LIB_DIR}
)