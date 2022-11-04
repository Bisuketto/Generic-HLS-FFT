# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(LUT_gen)

set(CMAKE_CXX_STANDARD 17)

add_executable(LUT_gen LUT_gen.cpp)

add_custom_command(
        TARGET LUT_gen
        POST_BUILD
        COMMAND LUT_gen ${CMAKE_CURRENT_SOURCE_DIR}/../../../roms
        COMMENT "COMPUTING ROMS"
        BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../roms/
        VERBATIM
)