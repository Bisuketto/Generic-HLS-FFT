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
 * @file dig_rev_perm_num_stage.h
 * @author Guanglin Xu (guanglix (at) xxxandrew.cmu.edu (delete xxx)), modified by Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains reformated digit reverse stage configuration from the SPIRAL project, originally written by Guanglin Xu
 * @version 0.0.0
 * @date 2022-01-17
 * 
 * @license This source is released under the SPIRAL License.
 * 
 */

#ifndef DIGIT_REV_PERM_NUM_STAGE_H_
#define DIGIT_REV_PERM_NUM_STAGE_H_

#define DIGIT_REV_NUM_STAGE_N16_SW2_R2      1
#define DIGIT_REV_NUM_STAGE_N16_SW4_R2      2
#define DIGIT_REV_NUM_STAGE_N16_SW8_R2      1
#define DIGIT_REV_NUM_STAGE_N32_SW2_R2      1
#define DIGIT_REV_NUM_STAGE_N32_SW4_R2      2
#define DIGIT_REV_NUM_STAGE_N32_SW8_R2      2
#define DIGIT_REV_NUM_STAGE_N64_SW2_R2      1
#define DIGIT_REV_NUM_STAGE_N64_SW4_R2      2
#define DIGIT_REV_NUM_STAGE_N64_SW8_R2      3
#define DIGIT_REV_NUM_STAGE_N128_SW2_R2     1
#define DIGIT_REV_NUM_STAGE_N128_SW4_R2     2
#define DIGIT_REV_NUM_STAGE_N128_SW8_R2     3
#define DIGIT_REV_NUM_STAGE_N256_SW2_R2     1
#define DIGIT_REV_NUM_STAGE_N256_SW4_R2     2
#define DIGIT_REV_NUM_STAGE_N256_SW8_R2     3
#define DIGIT_REV_NUM_STAGE_N512_SW2_R2     1
#define DIGIT_REV_NUM_STAGE_N512_SW4_R2     2
#define DIGIT_REV_NUM_STAGE_N512_SW8_R2     3
#define DIGIT_REV_NUM_STAGE_N1024_SW2_R2    1
#define DIGIT_REV_NUM_STAGE_N1024_SW4_R2    2
#define DIGIT_REV_NUM_STAGE_N1024_SW8_R2    3
#define DIGIT_REV_NUM_STAGE_N2048_SW2_R2    1
#define DIGIT_REV_NUM_STAGE_N2048_SW4_R2    2
#define DIGIT_REV_NUM_STAGE_N2048_SW8_R2    3
#define DIGIT_REV_NUM_STAGE_N4096_SW2_R2    1
#define DIGIT_REV_NUM_STAGE_N4096_SW4_R2    2
#define DIGIT_REV_NUM_STAGE_N4096_SW8_R2    3
#define DIGIT_REV_NUM_STAGE_N16_SW4_R4      2
#define DIGIT_REV_NUM_STAGE_N16_SW8_R4      1
#define DIGIT_REV_NUM_STAGE_N64_SW4_R4      2
#define DIGIT_REV_NUM_STAGE_N64_SW8_R4      2
#define DIGIT_REV_NUM_STAGE_N256_SW4_R4     2
#define DIGIT_REV_NUM_STAGE_N256_SW8_R4     3
#define DIGIT_REV_NUM_STAGE_N1024_SW4_R4    2
#define DIGIT_REV_NUM_STAGE_N1024_SW8_R4    3
#define DIGIT_REV_NUM_STAGE_N4096_SW4_R4    2
#define DIGIT_REV_NUM_STAGE_N4096_SW8_R4    3
#define DIGIT_REV_NUM_STAGE_N64_SW8_R8      3
#define DIGIT_REV_NUM_STAGE_N512_SW8_R8     3
#define DIGIT_REV_NUM_STAGE_N4096_SW8_R8    3


constexpr unsigned get_digit_rev_num_stage(unsigned N, unsigned SW, unsigned R){
    if( N == 16 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 16 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 16 && SW == 8 && R == 2){
        return 1;
    }
    else if( N == 32 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 32 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 32 && SW == 8 && R == 2){
        return 2;
    }
    else if( N == 64 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 64 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 64 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 128 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 128 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 128 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 256 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 256 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 256 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 512 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 512 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 512 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 1024 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 1024 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 1024 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 2048 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 2048 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 2048 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 4096 && SW == 2 && R == 2){
        return 1;
    }
    else if( N == 4096 && SW == 4 && R == 2){
        return 2;
    }
    else if( N == 4096 && SW == 8 && R == 2){
        return 3;
    }
    else if( N == 16 && SW == 4 && R == 4){
        return 2;
    }
    else if( N == 16 && SW == 8 && R == 4){
        return 1;
    }
    else if( N == 64 && SW == 4 && R == 4){
        return 2;
    }
    else if( N == 64 && SW == 8 && R == 4){
        return 2;
    }
    else if( N == 256 && SW == 4 && R == 4){
        return 2;
    }
    else if( N == 256 && SW == 8 && R == 4){
        return 3;
    }
    else if( N == 1024 && SW == 4 && R == 4){
        return 2;
    }
    else if( N == 1024 && SW == 8 && R == 4){
        return 3;
    }
    else if( N == 4096 && SW == 4 && R == 4){
        return 2;
    }
    else if( N == 4096 && SW == 8 && R == 4){
        return 3;
    }
    else if( N == 64 && SW == 8 && R == 8){
        return 3;
    }
    else if( N == 512 && SW == 8 && R == 8){
        return 3;
    }
    else if( N == 4096 && SW == 8 && R == 8){
        return 3;
    }
    else{
        return 0;
    }
}

#endif