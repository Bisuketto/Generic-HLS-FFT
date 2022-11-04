# Author : Hugues ALMORIN (hugues.almorin@arelis.com)
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3

cmake_minimum_required(VERSION 3.16.3)

project(SIG_gen)

set(CMAKE_CXX_STANDARD 17)

add_executable(SIG_gen SIG_gen.cpp)

add_custom_command(
        TARGET SIG_gen
        POST_BUILD
        COMMAND SIG_gen ${CMAKE_CURRENT_SOURCE_DIR}/../../../signals
        COMMENT "COMPUTING TEST SIGNALS"
        BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../signals/
        VERBATIM
)