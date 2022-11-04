// SPIRAL License
//
// Copyright 2017, Carnegie Mellon University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the SPIRAL project.

/**
 * @file stride_perm_config.h
 * @author Guanglin Xu (guanglix (at) xxxandrew.cmu.edu (delete xxx)), modified by Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains reformated stride permutation configuration from the SPIRAL project, originally written by Guanglin Xu
 * @version 0.0.0
 * @date 2022-01-17
 * 
 * @license This source is released under the SPIRAL License.
 * 
 */

#ifndef STRIDE_PERM_CONFIG_H_
#define STRIDE_PERM_CONFIG_H_

#define STRIDE_PERM_CONFIG_N16_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 3, 2}, {{0, 1}}, {3}
#define STRIDE_PERM_CONFIG_N16_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 3}, {{0, 2, 1, 3}}, {3}
#define STRIDE_PERM_CONFIG_N16_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}
#define STRIDE_PERM_CONFIG_N32_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 4, 3, 2}, {{0, 1}}, {4}
#define STRIDE_PERM_CONFIG_N32_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 4, 3}, {{0, 2, 1, 3}}, {4}
#define STRIDE_PERM_CONFIG_N32_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {4}
#define STRIDE_PERM_CONFIG_N64_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 5, 4, 3, 2}, {{0, 1}}, {5}
#define STRIDE_PERM_CONFIG_N64_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 5, 4, 3}, {{0, 2, 1, 3}}, {5}
#define STRIDE_PERM_CONFIG_N64_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {5}
#define STRIDE_PERM_CONFIG_N128_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 6, 5, 4, 3, 2}, {{0, 1}}, {6}
#define STRIDE_PERM_CONFIG_N128_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {6}
#define STRIDE_PERM_CONFIG_N128_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {6}
#define STRIDE_PERM_CONFIG_N256_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 7, 6, 5, 4, 3, 2}, {{0, 1}}, {7}
#define STRIDE_PERM_CONFIG_N256_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 7, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {7}
#define STRIDE_PERM_CONFIG_N256_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 7, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {7}
#define STRIDE_PERM_CONFIG_N512_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 8, 7, 6, 5, 4, 3, 2}, {{0, 1}}, {8}
#define STRIDE_PERM_CONFIG_N512_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 8, 7, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {8}
#define STRIDE_PERM_CONFIG_N512_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 8, 7, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {8}
#define STRIDE_PERM_CONFIG_N1024_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 9, 8, 7, 6, 5, 4, 3, 2}, {{0, 1}}, {9}
#define STRIDE_PERM_CONFIG_N1024_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 9, 8, 7, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {9}
#define STRIDE_PERM_CONFIG_N1024_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 9, 8, 7, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {9}
#define STRIDE_PERM_CONFIG_N2048_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 10, 9, 8, 7, 6, 5, 4, 3, 2}, {{0, 1}}, {10}
#define STRIDE_PERM_CONFIG_N2048_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 10, 9, 8, 7, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {10}
#define STRIDE_PERM_CONFIG_N2048_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 10, 9, 8, 7, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {10}
#define STRIDE_PERM_CONFIG_N4096_SW2_R2	    {0, 1}, {{0, 1}}, {1}, {0, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2}, {{0, 1}}, {11}
#define STRIDE_PERM_CONFIG_N4096_SW4_R2	    {0, 2, 1, 3}, {{0, 2, 1, 3}}, {2}, {0, 11, 10, 9, 8, 7, 6, 5, 4, 3}, {{0, 2, 1, 3}}, {11}
#define STRIDE_PERM_CONFIG_N4096_SW8_R2	    {0, 4, 1, 5, 2, 6, 3, 7}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {3}, {0, 11, 10, 9, 8, 7, 6, 5, 4}, {{0, 4, 1, 5, 2, 6, 3, 7}}, {11}
#define STRIDE_PERM_CONFIG_N16_SW4_R4	    {0, 1, 2, 3}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}, {1, 0}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}
#define STRIDE_PERM_CONFIG_N16_SW8_R4	    {0, 4, 2, 6, 1, 5, 3, 7}, {{0, 2, 1, 3, 4, 6, 5, 7}}, {3}, {1}, {{0, 2, 1, 3, 4, 6, 5, 7}}, {3}
#define STRIDE_PERM_CONFIG_N64_SW4_R4	    {0, 1, 2, 3}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}, {1, 0, 5, 4}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {4, 5}
#define STRIDE_PERM_CONFIG_N64_SW8_R4	    {0, 2, 4, 6, 1, 3, 5, 7}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4}, {1, 0, 5}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {4, 5}
#define STRIDE_PERM_CONFIG_N256_SW4_R4	    {0, 1, 2, 3}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}, {1, 0, 7, 6, 5, 4}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {6, 7}
#define STRIDE_PERM_CONFIG_N256_SW8_R4	    {0, 2, 4, 6, 1, 3, 5, 7}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4}, {1, 0, 7, 6, 5}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {6, 7}
#define STRIDE_PERM_CONFIG_N1024_SW4_R4	    {0, 1, 2, 3}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}, {1, 0, 9, 8, 7, 6, 5, 4}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {8, 9}
#define STRIDE_PERM_CONFIG_N1024_SW8_R4	    {0, 2, 4, 6, 1, 3, 5, 7}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4}, {1, 0, 9, 8, 7, 6, 5}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {8, 9}
#define STRIDE_PERM_CONFIG_N4096_SW4_R4	    {0, 1, 2, 3}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {2, 3}, {1, 0, 11, 10, 9, 8, 7, 6, 5, 4}, {{0, 1, 2, 3},{0, 2, 1, 3}}, {10, 11}
#define STRIDE_PERM_CONFIG_N4096_SW8_R4	    {0, 2, 4, 6, 1, 3, 5, 7}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4}, {1, 0, 11, 10, 9, 8, 7, 6, 5}, {{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {10, 11}
#define STRIDE_PERM_CONFIG_N64_SW8_R8	    {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4, 5}, {2, 1, 0}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4, 5}
#define STRIDE_PERM_CONFIG_N512_SW8_R8	    {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4, 5}, {2, 1, 0, 8, 7, 6}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {6, 7, 8}
#define STRIDE_PERM_CONFIG_N4096_SW8_R8	    {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {3, 4, 5}, {2, 1, 0, 11, 10, 9, 8, 7, 6}, {{0, 1, 2, 3, 4, 5, 6, 7},{0, 2, 1, 3, 4, 6, 5, 7},{0, 4, 1, 5, 2, 6, 3, 7}}, {9, 10, 11}

#endif