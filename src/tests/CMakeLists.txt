# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(FFT_calc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(FFT_calc fft_calc.cpp)

add_dependencies(FFT_calc SIG_gen)

target_include_directories( FFT_calc PRIVATE 
        ${AP_TYPES_LIB_DIR}
)

target_link_libraries(FFT_calc fftw3
        fft_fftw
        fft_sbg_radix
        fft_tools
)

add_custom_command(
        TARGET FFT_calc
        POST_BUILD
        COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/../../data/results
        COMMENT "CREATING results FOLDER"
        BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/../../data/results/ ${CMAKE_CURRENT_SOURCE_DIR}/../../data/
        VERBATIM
)