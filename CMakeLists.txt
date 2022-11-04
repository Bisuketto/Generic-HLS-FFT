# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)
include(CMakePrintHelpers)
set(AP_TYPES_LIB_DIR
    /mnt/c/Users/halmorin/Documents/THESE/dev/HLS_arbitrary_Precision_Types/include
    CACHE PATH
    "Path to ap_types library directory for simulation"
)
cmake_print_variables(AP_TYPES_LIB_DIR)

project(Generic-HLS-FFT LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src/common)
add_subdirectory(src/fft)
add_subdirectory(src/fft_fftw)
add_subdirectory(src/tests)
add_subdirectory(src/tools/LUT_gen)
add_subdirectory(src/tools/SIG_gen)
