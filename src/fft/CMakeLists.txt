# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(fft_sbg_radix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(fft_sbg_radix sw_wrap/lib_fft_sbg_radix.cpp)
target_include_directories( fft_sbg_radix PRIVATE 
        ${AP_TYPES_LIB_DIR}
)
target_link_libraries(fft_sbg_radix fft_tools)

add_dependencies(fft_sbg_radix LUT_gen)