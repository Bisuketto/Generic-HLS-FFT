/**
 * @file lib_fft_sbg_radix.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains software wrapping for using the FFT model
 * @version 0.0.0
 * @date 2022-03-10
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include "lib_fft_sbg_radix.hpp"

static bool rom_loaded = false;
static float* t_cos;
static float* t_sin;
static int64_t* t_cos_fixed;
static int64_t* t_sin_fixed;

void lib_fft_sbg_radix_load_fft_roms(std::string cfile, const int32_t SIZE){
    t_cos = new float[SIZE];
    t_sin = new float[SIZE];
    t_cos_fixed = new int64_t[SIZE];
    t_sin_fixed = new int64_t[SIZE];

    std::string file_cos = cfile + "/cos.txt";
    std::string file_sin = cfile + "/sin.txt";

    //
    // CHARGEMENT DE LA TABLE DES COSINUS
    //

    std::ifstream cosfile( file_cos );
    if( cosfile.is_open() == false ){
        printf("(EE) Probleme à l'ouverture du fichier (%s)\n", file_cos.c_str());
        exit( EXIT_FAILURE );
    }
    for(int i = 0; i < SIZE; i += 1)
    {
        std::string tmp;
        std::getline(cosfile, tmp);
        t_cos[i] = std::stof( tmp );
        t_cos_fixed[i] = std::stoi( tmp );
    }
    cosfile.close();

    //
    // CHARGEMENT DE LA TABLE DES SINUS
    //

    std::ifstream sinfile( file_sin );
    if( sinfile.is_open() == false ){
        printf("(EE) Probleme à l'ouverture du fichier (%s)\n", file_sin.c_str());
        exit( EXIT_FAILURE );
    }
    for(int i = 0; i < SIZE; i += 1)
    {
        std::string tmp;
        std::getline(sinfile, tmp);
        t_sin[i] = std::stof( tmp );
        t_sin_fixed[i] = std::stoi( tmp );
        
        t_sin[i] = t_sin[i]; // On inverse la table du sinus car dans l'algo on va de -2pi à 0
                              // Et dans le fichier on fait l'inverse !
        t_sin_fixed[i] = t_sin_fixed[i];
    }
    sinfile.close();

    rom_loaded = true;
}

void lib_sbg_radix_fft(DataVector& data, Parameters& p){
    const std::string q_input    = p.toString("q_input");
    const std::string q_internal = p.toString("q_internal");
    const std::string q_rom      = p.toString("q_rom");

    const int32_t SIZE = p.toInt("FFT_SIZE");
    const int32_t SW = p.toInt("SW");
    const int32_t RADIX = p.toInt("RADIX");
    // const int32_t RADIX = 4;
    // const int32_t SW = 8;
    constexpr const int GS = 2;

    if(!rom_loaded){
        std::cout << "(EE) Kastner Lut ROM not loaded !" << std::endl;
        exit( EXIT_FAILURE );
    }

    if( (q_internal=="float") && (q_rom=="float") )
    {   
        float* inI = new float[SIZE];
        float* inQ = new float[SIZE];
        float* outI = new float[SIZE];
        float* outQ = new float[SIZE];

        if((q_input == "double") || (q_input == "float")){
            for(int i = 0; i < SIZE; i++){
                inI[i] = data.I[i];
                inQ[i] = data.Q[i];//samples_in_Q.read();
            }
        }
        else{
            const uint8_t Q_in = p.toInt("q_input");
            for(int i = 0; i < SIZE; i++){
                inI[i] = std::round(data.I[i]*(1 << (Q_in-1)));
                inQ[i] = std::round(data.Q[i]*(1 << (Q_in-1)));//samples_in_Q.read();
            }
        }

        perm_config<DIGIT_REV_NUM_STAGE_N16_SW2_R2  , 2, ceillog2(  16), ceillog2(2)> dig_rev_config_R2_N16_SW2   = { DIG_REV_PERM_CONFIG_N16_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW2_R2  , 2, ceillog2(  32), ceillog2(2)> dig_rev_config_R2_N32_SW2   = { DIG_REV_PERM_CONFIG_N32_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW2_R2  , 2, ceillog2(  64), ceillog2(2)> dig_rev_config_R2_N64_SW2   = { DIG_REV_PERM_CONFIG_N64_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW2_R2 , 2, ceillog2( 128), ceillog2(2)> dig_rev_config_R2_N128_SW2  = { DIG_REV_PERM_CONFIG_N128_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW2_R2 , 2, ceillog2( 256), ceillog2(2)> dig_rev_config_R2_N256_SW2  = { DIG_REV_PERM_CONFIG_N256_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW2_R2 , 2, ceillog2( 512), ceillog2(2)> dig_rev_config_R2_N512_SW2  = { DIG_REV_PERM_CONFIG_N512_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW2_R2, 2, ceillog2(1024), ceillog2(2)> dig_rev_config_R2_N1024_SW2 = { DIG_REV_PERM_CONFIG_N1024_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW2_R2, 2, ceillog2(2048), ceillog2(2)> dig_rev_config_R2_N2048_SW2 = { DIG_REV_PERM_CONFIG_N2048_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW2_R2, 2, ceillog2(4096), ceillog2(2)> dig_rev_config_R2_N4096_SW2 = { DIG_REV_PERM_CONFIG_N4096_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW4_R2  , 4, ceillog2(  16), ceillog2(4)> dig_rev_config_R2_N16_SW4   = { DIG_REV_PERM_CONFIG_N16_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW4_R2  , 4, ceillog2(  32), ceillog2(4)> dig_rev_config_R2_N32_SW4   = { DIG_REV_PERM_CONFIG_N32_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW4_R2  , 4, ceillog2(  64), ceillog2(4)> dig_rev_config_R2_N64_SW4   = { DIG_REV_PERM_CONFIG_N64_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW4_R2 , 4, ceillog2( 128), ceillog2(4)> dig_rev_config_R2_N128_SW4  = { DIG_REV_PERM_CONFIG_N128_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW4_R2 , 4, ceillog2( 256), ceillog2(4)> dig_rev_config_R2_N256_SW4  = { DIG_REV_PERM_CONFIG_N256_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW4_R2 , 4, ceillog2( 512), ceillog2(4)> dig_rev_config_R2_N512_SW4  = { DIG_REV_PERM_CONFIG_N512_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW4_R2, 4, ceillog2(1024), ceillog2(4)> dig_rev_config_R2_N1024_SW4 = { DIG_REV_PERM_CONFIG_N1024_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW4_R2, 4, ceillog2(2048), ceillog2(4)> dig_rev_config_R2_N2048_SW4 = { DIG_REV_PERM_CONFIG_N2048_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW4_R2, 4, ceillog2(4096), ceillog2(4)> dig_rev_config_R2_N4096_SW4 = { DIG_REV_PERM_CONFIG_N4096_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW8_R2  , 8, ceillog2(  16), ceillog2(8)> dig_rev_config_R2_N16_SW8   = { DIG_REV_PERM_CONFIG_N16_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW8_R2  , 8, ceillog2(  32), ceillog2(8)> dig_rev_config_R2_N32_SW8   = { DIG_REV_PERM_CONFIG_N32_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R2  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R2_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW8_R2 , 8, ceillog2( 128), ceillog2(8)> dig_rev_config_R2_N128_SW8  = { DIG_REV_PERM_CONFIG_N128_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW8_R2 , 8, ceillog2( 256), ceillog2(8)> dig_rev_config_R2_N256_SW8  = { DIG_REV_PERM_CONFIG_N256_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW8_R2 , 8, ceillog2( 512), ceillog2(8)> dig_rev_config_R2_N512_SW8  = { DIG_REV_PERM_CONFIG_N512_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW8_R2, 8, ceillog2(1024), ceillog2(8)> dig_rev_config_R2_N1024_SW8 = { DIG_REV_PERM_CONFIG_N1024_SW8_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW8_R2, 8, ceillog2(2048), ceillog2(8)> dig_rev_config_R2_N2048_SW8 = { DIG_REV_PERM_CONFIG_N2048_SW8_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R2, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R2_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R2     };

        perm_config<DIGIT_REV_NUM_STAGE_N16_SW4_R4  , 4, ceillog2(  16), ceillog2(4)> dig_rev_config_R4_N16_SW4   = { DIG_REV_PERM_CONFIG_N16_SW4_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW4_R4  , 4, ceillog2(  64), ceillog2(4)> dig_rev_config_R4_N64_SW4   = { DIG_REV_PERM_CONFIG_N64_SW4_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW4_R4 , 4, ceillog2( 256), ceillog2(4)> dig_rev_config_R4_N256_SW4  = { DIG_REV_PERM_CONFIG_N256_SW4_R4      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW4_R4, 4, ceillog2(1024), ceillog2(4)> dig_rev_config_R4_N1024_SW4 = { DIG_REV_PERM_CONFIG_N1024_SW4_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW4_R4, 4, ceillog2(4096), ceillog2(4)> dig_rev_config_R4_N4096_SW4 = { DIG_REV_PERM_CONFIG_N4096_SW4_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW8_R4  , 8, ceillog2(  16), ceillog2(8)> dig_rev_config_R4_N16_SW8   = { DIG_REV_PERM_CONFIG_N16_SW8_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R4  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R4_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW8_R4 , 8, ceillog2( 256), ceillog2(8)> dig_rev_config_R4_N256_SW8  = { DIG_REV_PERM_CONFIG_N256_SW8_R4      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW8_R4, 8, ceillog2(1024), ceillog2(8)> dig_rev_config_R4_N1024_SW8 = { DIG_REV_PERM_CONFIG_N1024_SW8_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R4, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R4_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R4     };

        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R8  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R8_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R8       };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW8_R8 , 8, ceillog2( 512), ceillog2(8)> dig_rev_config_R8_N512_SW8  = { DIG_REV_PERM_CONFIG_N512_SW8_R8      };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R8, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R8_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R8     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW2_R2  , 2, ceillog2(  16), ceillog2(2)> stride_config_R2_N16_SW2   = { STRIDE_PERM_CONFIG_N16_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW2_R2  , 2, ceillog2(  32), ceillog2(2)> stride_config_R2_N32_SW2   = { STRIDE_PERM_CONFIG_N32_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW2_R2  , 2, ceillog2(  64), ceillog2(2)> stride_config_R2_N64_SW2   = { STRIDE_PERM_CONFIG_N64_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW2_R2 , 2, ceillog2( 128), ceillog2(2)> stride_config_R2_N128_SW2  = { STRIDE_PERM_CONFIG_N128_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW2_R2 , 2, ceillog2( 256), ceillog2(2)> stride_config_R2_N256_SW2  = { STRIDE_PERM_CONFIG_N256_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW2_R2 , 2, ceillog2( 512), ceillog2(2)> stride_config_R2_N512_SW2  = { STRIDE_PERM_CONFIG_N512_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW2_R2, 2, ceillog2(1024), ceillog2(2)> stride_config_R2_N1024_SW2 = { STRIDE_PERM_CONFIG_N1024_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW2_R2, 2, ceillog2(2048), ceillog2(2)> stride_config_R2_N2048_SW2 = { STRIDE_PERM_CONFIG_N2048_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW2_R2, 2, ceillog2(4096), ceillog2(2)> stride_config_R2_N4096_SW2 = { STRIDE_PERM_CONFIG_N4096_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R2  , 4, ceillog2(  16), ceillog2(4)> stride_config_R2_N16_SW4   = { STRIDE_PERM_CONFIG_N16_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW4_R2  , 4, ceillog2(  32), ceillog2(4)> stride_config_R2_N32_SW4   = { STRIDE_PERM_CONFIG_N32_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R2  , 4, ceillog2(  64), ceillog2(4)> stride_config_R2_N64_SW4   = { STRIDE_PERM_CONFIG_N64_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW4_R2 , 4, ceillog2( 128), ceillog2(4)> stride_config_R2_N128_SW4  = { STRIDE_PERM_CONFIG_N128_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R2 , 4, ceillog2( 256), ceillog2(4)> stride_config_R2_N256_SW4  = { STRIDE_PERM_CONFIG_N256_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW4_R2 , 4, ceillog2( 512), ceillog2(4)> stride_config_R2_N512_SW4  = { STRIDE_PERM_CONFIG_N512_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R2, 4, ceillog2(1024), ceillog2(4)> stride_config_R2_N1024_SW4 = { STRIDE_PERM_CONFIG_N1024_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW4_R2, 4, ceillog2(2048), ceillog2(4)> stride_config_R2_N2048_SW4 = { STRIDE_PERM_CONFIG_N2048_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R2, 4, ceillog2(4096), ceillog2(4)> stride_config_R2_N4096_SW4 = { STRIDE_PERM_CONFIG_N4096_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R2  , 8, ceillog2(  16), ceillog2(8)> stride_config_R2_N16_SW8   = { STRIDE_PERM_CONFIG_N16_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW8_R2  , 8, ceillog2(  32), ceillog2(8)> stride_config_R2_N32_SW8   = { STRIDE_PERM_CONFIG_N32_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R2  , 8, ceillog2(  64), ceillog2(8)> stride_config_R2_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW8_R2 , 8, ceillog2( 128), ceillog2(8)> stride_config_R2_N128_SW8  = { STRIDE_PERM_CONFIG_N128_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R2 , 8, ceillog2( 256), ceillog2(8)> stride_config_R2_N256_SW8  = { STRIDE_PERM_CONFIG_N256_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R2 , 8, ceillog2( 512), ceillog2(8)> stride_config_R2_N512_SW8  = { STRIDE_PERM_CONFIG_N512_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R2, 8, ceillog2(1024), ceillog2(8)> stride_config_R2_N1024_SW8 = { STRIDE_PERM_CONFIG_N1024_SW8_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW8_R2, 8, ceillog2(2048), ceillog2(8)> stride_config_R2_N2048_SW8 = { STRIDE_PERM_CONFIG_N2048_SW8_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R2, 8, ceillog2(4096), ceillog2(8)> stride_config_R2_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R2     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R4  , 4, ceillog2(  16), ceillog2(4)> stride_config_R4_N16_SW4   = { STRIDE_PERM_CONFIG_N16_SW4_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R4  , 4, ceillog2(  64), ceillog2(4)> stride_config_R4_N64_SW4   = { STRIDE_PERM_CONFIG_N64_SW4_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R4 , 4, ceillog2( 256), ceillog2(4)> stride_config_R4_N256_SW4  = { STRIDE_PERM_CONFIG_N256_SW4_R4      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R4, 4, ceillog2(1024), ceillog2(4)> stride_config_R4_N1024_SW4 = { STRIDE_PERM_CONFIG_N1024_SW4_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R4, 4, ceillog2(4096), ceillog2(4)> stride_config_R4_N4096_SW4 = { STRIDE_PERM_CONFIG_N4096_SW4_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R4  , 8, ceillog2(  16), ceillog2(8)> stride_config_R4_N16_SW8   = { STRIDE_PERM_CONFIG_N16_SW8_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R4  , 8, ceillog2(  64), ceillog2(8)> stride_config_R4_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R4 , 8, ceillog2( 256), ceillog2(8)> stride_config_R4_N256_SW8  = { STRIDE_PERM_CONFIG_N256_SW8_R4      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R4, 8, ceillog2(1024), ceillog2(8)> stride_config_R4_N1024_SW8 = { STRIDE_PERM_CONFIG_N1024_SW8_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R4, 8, ceillog2(4096), ceillog2(8)> stride_config_R4_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R4     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R8  , 8, ceillog2(  64), ceillog2(8)> stride_config_R8_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R8       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R8 , 8, ceillog2( 512), ceillog2(8)> stride_config_R8_N512_SW8  = { STRIDE_PERM_CONFIG_N512_SW8_R8      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R8, 8, ceillog2(4096), ceillog2(8)> stride_config_R8_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R8     };
        
        if(RADIX == 2){
            if(SW == 2){
                switch (SIZE) {
                    case   16: sbg_radix_fft<float, float, float,   16, 2, 2, GS, ceillog2(  16), ceillog2(2), DIGIT_REV_NUM_STAGE_N16_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW2_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N16_SW2  , stride_config_R2_N16_SW2  , 0); break;
                    case   32: sbg_radix_fft<float, float, float,   32, 2, 2, GS, ceillog2(  32), ceillog2(2), DIGIT_REV_NUM_STAGE_N32_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW2_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N32_SW2  , stride_config_R2_N32_SW2  , 0); break;
                    case   64: sbg_radix_fft<float, float, float,   64, 2, 2, GS, ceillog2(  64), ceillog2(2), DIGIT_REV_NUM_STAGE_N64_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW2_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N64_SW2  , stride_config_R2_N64_SW2  , 0); break;
                    case  128: sbg_radix_fft<float, float, float,  128, 2, 2, GS, ceillog2( 128), ceillog2(2), DIGIT_REV_NUM_STAGE_N128_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW2_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N128_SW2 , stride_config_R2_N128_SW2 , 0); break;
                    case  256: sbg_radix_fft<float, float, float,  256, 2, 2, GS, ceillog2( 256), ceillog2(2), DIGIT_REV_NUM_STAGE_N256_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW2_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N256_SW2 , stride_config_R2_N256_SW2 , 0); break;
                    case  512: sbg_radix_fft<float, float, float,  512, 2, 2, GS, ceillog2( 512), ceillog2(2), DIGIT_REV_NUM_STAGE_N512_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW2_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N512_SW2 , stride_config_R2_N512_SW2 , 0); break;
                    case 1024: sbg_radix_fft<float, float, float, 1024, 2, 2, GS, ceillog2(1024), ceillog2(2), DIGIT_REV_NUM_STAGE_N1024_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW2_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N1024_SW2, stride_config_R2_N1024_SW2, 0); break;
                    case 2048: sbg_radix_fft<float, float, float, 2048, 2, 2, GS, ceillog2(2048), ceillog2(2), DIGIT_REV_NUM_STAGE_N2048_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW2_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N2048_SW2, stride_config_R2_N2048_SW2, 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 2, 2, GS, ceillog2(4096), ceillog2(2), DIGIT_REV_NUM_STAGE_N4096_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW2_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N4096_SW2, stride_config_R2_N4096_SW2, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }else if(SW == 4){
                switch (SIZE) {
                    case   16: sbg_radix_fft<float, float, float,   16, 2, 4, GS, ceillog2(  16), ceillog2(4), DIGIT_REV_NUM_STAGE_N16_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N16_SW4  , stride_config_R2_N16_SW4  , 0); break;
                    case   32: sbg_radix_fft<float, float, float,   32, 2, 4, GS, ceillog2(  32), ceillog2(4), DIGIT_REV_NUM_STAGE_N32_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW4_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N32_SW4  , stride_config_R2_N32_SW4  , 0); break;
                    case   64: sbg_radix_fft<float, float, float,   64, 2, 4, GS, ceillog2(  64), ceillog2(4), DIGIT_REV_NUM_STAGE_N64_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N64_SW4  , stride_config_R2_N64_SW4  , 0); break;
                    case  128: sbg_radix_fft<float, float, float,  128, 2, 4, GS, ceillog2( 128), ceillog2(4), DIGIT_REV_NUM_STAGE_N128_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW4_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N128_SW4 , stride_config_R2_N128_SW4 , 0); break;
                    case  256: sbg_radix_fft<float, float, float,  256, 2, 4, GS, ceillog2( 256), ceillog2(4), DIGIT_REV_NUM_STAGE_N256_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N256_SW4 , stride_config_R2_N256_SW4 , 0); break;
                    case  512: sbg_radix_fft<float, float, float,  512, 2, 4, GS, ceillog2( 512), ceillog2(4), DIGIT_REV_NUM_STAGE_N512_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW4_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N512_SW4 , stride_config_R2_N512_SW4 , 0); break;
                    case 1024: sbg_radix_fft<float, float, float, 1024, 2, 4, GS, ceillog2(1024), ceillog2(4), DIGIT_REV_NUM_STAGE_N1024_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N1024_SW4, stride_config_R2_N1024_SW4, 0); break;
                    case 2048: sbg_radix_fft<float, float, float, 2048, 2, 4, GS, ceillog2(2048), ceillog2(4), DIGIT_REV_NUM_STAGE_N2048_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW4_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N2048_SW4, stride_config_R2_N2048_SW4, 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 2, 4, GS, ceillog2(4096), ceillog2(4), DIGIT_REV_NUM_STAGE_N4096_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N4096_SW4, stride_config_R2_N4096_SW4, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            } else if(SW == 8){
                switch (SIZE) {
                    case   16: sbg_radix_fft<float, float, float,   16, 2, 8, GS, ceillog2(  16), ceillog2(8), DIGIT_REV_NUM_STAGE_N16_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N16_SW8  , stride_config_R2_N16_SW8  , 0); break;
                    case   32: sbg_radix_fft<float, float, float,   32, 2, 8, GS, ceillog2(  32), ceillog2(8), DIGIT_REV_NUM_STAGE_N32_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW8_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N32_SW8  , stride_config_R2_N32_SW8  , 0); break;
                    case   64: sbg_radix_fft<float, float, float,   64, 2, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R2  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N64_SW8  , stride_config_R2_N64_SW8  , 0); break;
                    case  128: sbg_radix_fft<float, float, float,  128, 2, 8, GS, ceillog2( 128), ceillog2(8), DIGIT_REV_NUM_STAGE_N128_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW8_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N128_SW8 , stride_config_R2_N128_SW8 , 0); break;
                    case  256: sbg_radix_fft<float, float, float,  256, 2, 8, GS, ceillog2( 256), ceillog2(8), DIGIT_REV_NUM_STAGE_N256_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N256_SW8 , stride_config_R2_N256_SW8 , 0); break;
                    case  512: sbg_radix_fft<float, float, float,  512, 2, 8, GS, ceillog2( 512), ceillog2(8), DIGIT_REV_NUM_STAGE_N512_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R2 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N512_SW8 , stride_config_R2_N512_SW8 , 0); break;
                    case 1024: sbg_radix_fft<float, float, float, 1024, 2, 8, GS, ceillog2(1024), ceillog2(8), DIGIT_REV_NUM_STAGE_N1024_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N1024_SW8, stride_config_R2_N1024_SW8, 0); break;
                    case 2048: sbg_radix_fft<float, float, float, 2048, 2, 8, GS, ceillog2(2048), ceillog2(8), DIGIT_REV_NUM_STAGE_N2048_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW8_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N2048_SW8, stride_config_R2_N2048_SW8, 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 2, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R2>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R2_N4096_SW8, stride_config_R2_N4096_SW8, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }else if(RADIX == 4){
            if(SW == 4){
                switch (SIZE) {
                    case   16: sbg_radix_fft<float, float, float,   16, 4, 4, GS, ceillog2(  16), ceillog2(4), DIGIT_REV_NUM_STAGE_N16_SW4_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R4  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N16_SW4  , stride_config_R4_N16_SW4  , 0); break;
                    case   64: sbg_radix_fft<float, float, float,   64, 4, 4, GS, ceillog2(  64), ceillog2(4), DIGIT_REV_NUM_STAGE_N64_SW4_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R4  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N64_SW4  , stride_config_R4_N64_SW4  , 0); break;
                    case  256: sbg_radix_fft<float, float, float,  256, 4, 4, GS, ceillog2( 256), ceillog2(4), DIGIT_REV_NUM_STAGE_N256_SW4_R4 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R4 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N256_SW4 , stride_config_R4_N256_SW4 , 0); break;
                    case 1024: sbg_radix_fft<float, float, float, 1024, 4, 4, GS, ceillog2(1024), ceillog2(4), DIGIT_REV_NUM_STAGE_N1024_SW4_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R4>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N1024_SW4, stride_config_R4_N1024_SW4, 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 4, 4, GS, ceillog2(4096), ceillog2(4), DIGIT_REV_NUM_STAGE_N4096_SW4_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R4>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N4096_SW4, stride_config_R4_N4096_SW4, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            } else if(SW == 8){
                switch (SIZE) {
                    case   16: sbg_radix_fft<float, float, float,   16, 4, 8, GS, ceillog2(  16), ceillog2(8), DIGIT_REV_NUM_STAGE_N16_SW8_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R4  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N16_SW8  , stride_config_R4_N16_SW8  , 0); break;
                    case   64: sbg_radix_fft<float, float, float,   64, 4, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R4  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N64_SW8  , stride_config_R4_N64_SW8  , 0); break;
                    case  256: sbg_radix_fft<float, float, float,  256, 4, 8, GS, ceillog2( 256), ceillog2(8), DIGIT_REV_NUM_STAGE_N256_SW8_R4 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R4 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N256_SW8 , stride_config_R4_N256_SW8 , 0); break;
                    case 1024: sbg_radix_fft<float, float, float, 1024, 4, 8, GS, ceillog2(1024), ceillog2(8), DIGIT_REV_NUM_STAGE_N1024_SW8_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R4>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N1024_SW8, stride_config_R4_N1024_SW8, 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 4, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R4>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R4_N4096_SW8, stride_config_R4_N4096_SW8, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }else if(RADIX == 8){
            if(SW == 8){
                switch (SIZE) {
                    case   64: sbg_radix_fft<float, float, float,   64, 8, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R8  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R8  >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R8_N64_SW8  , stride_config_R8_N64_SW8  , 0); break;
                    case  512: sbg_radix_fft<float, float, float,  512, 8, 8, GS, ceillog2( 512), ceillog2(8), DIGIT_REV_NUM_STAGE_N512_SW8_R8 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R8 >(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R8_N512_SW8 , stride_config_R8_N512_SW8 , 0); break;
                    case 4096: sbg_radix_fft<float, float, float, 4096, 8, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R8, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R8>(inI, inQ, outI, outQ, t_cos, t_sin, dig_rev_config_R8_N4096_SW8, stride_config_R8_N4096_SW8, 0); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }
        
        
        if((q_input == "double") || (q_input == "float")){
            for(int i = 0; i < SIZE; i++){
                data.I[i]=outI[i];
                data.Q[i]=outQ[i];
            }
        }
        else{    
            const uint8_t Q_in = p.toInt("q_input");
            for(int i = 0; i < SIZE; i++){
                data.I[i]=outI[i]/(1<<(Q_in-1));
                data.Q[i]=outQ[i]/(1<<(Q_in-1));
            }
        }

        delete[] inI;
        delete[] inQ;
        delete[] outI;
        delete[] outQ;
    }else if((!((q_input == "double") || (q_input == "float"))) && (!((q_rom == "float") || (q_rom == "double"))))
    {
        const uint8_t Q_twiddle = p.toInt("q_rom"); // Q = QN.D with N = 2
        const uint8_t Q_in = p.toInt("q_input")-2;
        constexpr const int GS = 2;
        
        typedef ap_int<32> dType;
        typedef ap_int<32> tType;
        typedef ap_int<64> iType;
        
        dType* dinI = new dType[SIZE];
        dType* dinQ = new dType[SIZE];

        dType* outI = new dType[SIZE];
        dType* outQ = new dType[SIZE];


        tType* tc = new tType[SIZE];
        tType* ts = new tType[SIZE];

        for(int i = 0; i < SIZE; i++){
            dinI[i] = data.I.data()[i] * (((dType) 1) << Q_in);
            dinQ[i] = data.Q.data()[i] * (((dType) 1) << Q_in);
        }

        for(int i = 0; i < SIZE; i++){
            tc[i] = t_cos_fixed[i];
            ts[i] = t_sin_fixed[i];
        }

        perm_config<DIGIT_REV_NUM_STAGE_N16_SW2_R2  , 2, ceillog2(  16), ceillog2(2)> dig_rev_config_R2_N16_SW2   = { DIG_REV_PERM_CONFIG_N16_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW2_R2  , 2, ceillog2(  32), ceillog2(2)> dig_rev_config_R2_N32_SW2   = { DIG_REV_PERM_CONFIG_N32_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW2_R2  , 2, ceillog2(  64), ceillog2(2)> dig_rev_config_R2_N64_SW2   = { DIG_REV_PERM_CONFIG_N64_SW2_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW2_R2 , 2, ceillog2( 128), ceillog2(2)> dig_rev_config_R2_N128_SW2  = { DIG_REV_PERM_CONFIG_N128_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW2_R2 , 2, ceillog2( 256), ceillog2(2)> dig_rev_config_R2_N256_SW2  = { DIG_REV_PERM_CONFIG_N256_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW2_R2 , 2, ceillog2( 512), ceillog2(2)> dig_rev_config_R2_N512_SW2  = { DIG_REV_PERM_CONFIG_N512_SW2_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW2_R2, 2, ceillog2(1024), ceillog2(2)> dig_rev_config_R2_N1024_SW2 = { DIG_REV_PERM_CONFIG_N1024_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW2_R2, 2, ceillog2(2048), ceillog2(2)> dig_rev_config_R2_N2048_SW2 = { DIG_REV_PERM_CONFIG_N2048_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW2_R2, 2, ceillog2(4096), ceillog2(2)> dig_rev_config_R2_N4096_SW2 = { DIG_REV_PERM_CONFIG_N4096_SW2_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW4_R2  , 4, ceillog2(  16), ceillog2(4)> dig_rev_config_R2_N16_SW4   = { DIG_REV_PERM_CONFIG_N16_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW4_R2  , 4, ceillog2(  32), ceillog2(4)> dig_rev_config_R2_N32_SW4   = { DIG_REV_PERM_CONFIG_N32_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW4_R2  , 4, ceillog2(  64), ceillog2(4)> dig_rev_config_R2_N64_SW4   = { DIG_REV_PERM_CONFIG_N64_SW4_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW4_R2 , 4, ceillog2( 128), ceillog2(4)> dig_rev_config_R2_N128_SW4  = { DIG_REV_PERM_CONFIG_N128_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW4_R2 , 4, ceillog2( 256), ceillog2(4)> dig_rev_config_R2_N256_SW4  = { DIG_REV_PERM_CONFIG_N256_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW4_R2 , 4, ceillog2( 512), ceillog2(4)> dig_rev_config_R2_N512_SW4  = { DIG_REV_PERM_CONFIG_N512_SW4_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW4_R2, 4, ceillog2(1024), ceillog2(4)> dig_rev_config_R2_N1024_SW4 = { DIG_REV_PERM_CONFIG_N1024_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW4_R2, 4, ceillog2(2048), ceillog2(4)> dig_rev_config_R2_N2048_SW4 = { DIG_REV_PERM_CONFIG_N2048_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW4_R2, 4, ceillog2(4096), ceillog2(4)> dig_rev_config_R2_N4096_SW4 = { DIG_REV_PERM_CONFIG_N4096_SW4_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW8_R2  , 8, ceillog2(  16), ceillog2(8)> dig_rev_config_R2_N16_SW8   = { DIG_REV_PERM_CONFIG_N16_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N32_SW8_R2  , 8, ceillog2(  32), ceillog2(8)> dig_rev_config_R2_N32_SW8   = { DIG_REV_PERM_CONFIG_N32_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R2  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R2_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R2       };
        perm_config<DIGIT_REV_NUM_STAGE_N128_SW8_R2 , 8, ceillog2( 128), ceillog2(8)> dig_rev_config_R2_N128_SW8  = { DIG_REV_PERM_CONFIG_N128_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW8_R2 , 8, ceillog2( 256), ceillog2(8)> dig_rev_config_R2_N256_SW8  = { DIG_REV_PERM_CONFIG_N256_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW8_R2 , 8, ceillog2( 512), ceillog2(8)> dig_rev_config_R2_N512_SW8  = { DIG_REV_PERM_CONFIG_N512_SW8_R2      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW8_R2, 8, ceillog2(1024), ceillog2(8)> dig_rev_config_R2_N1024_SW8 = { DIG_REV_PERM_CONFIG_N1024_SW8_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N2048_SW8_R2, 8, ceillog2(2048), ceillog2(8)> dig_rev_config_R2_N2048_SW8 = { DIG_REV_PERM_CONFIG_N2048_SW8_R2     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R2, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R2_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R2     };

        perm_config<DIGIT_REV_NUM_STAGE_N16_SW4_R4  , 4, ceillog2(  16), ceillog2(4)> dig_rev_config_R4_N16_SW4   = { DIG_REV_PERM_CONFIG_N16_SW4_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW4_R4  , 4, ceillog2(  64), ceillog2(4)> dig_rev_config_R4_N64_SW4   = { DIG_REV_PERM_CONFIG_N64_SW4_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW4_R4 , 4, ceillog2( 256), ceillog2(4)> dig_rev_config_R4_N256_SW4  = { DIG_REV_PERM_CONFIG_N256_SW4_R4      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW4_R4, 4, ceillog2(1024), ceillog2(4)> dig_rev_config_R4_N1024_SW4 = { DIG_REV_PERM_CONFIG_N1024_SW4_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW4_R4, 4, ceillog2(4096), ceillog2(4)> dig_rev_config_R4_N4096_SW4 = { DIG_REV_PERM_CONFIG_N4096_SW4_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N16_SW8_R4  , 8, ceillog2(  16), ceillog2(8)> dig_rev_config_R4_N16_SW8   = { DIG_REV_PERM_CONFIG_N16_SW8_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R4  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R4_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R4       };
        perm_config<DIGIT_REV_NUM_STAGE_N256_SW8_R4 , 8, ceillog2( 256), ceillog2(8)> dig_rev_config_R4_N256_SW8  = { DIG_REV_PERM_CONFIG_N256_SW8_R4      };
        perm_config<DIGIT_REV_NUM_STAGE_N1024_SW8_R4, 8, ceillog2(1024), ceillog2(8)> dig_rev_config_R4_N1024_SW8 = { DIG_REV_PERM_CONFIG_N1024_SW8_R4     };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R4, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R4_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R4     };

        perm_config<DIGIT_REV_NUM_STAGE_N64_SW8_R8  , 8, ceillog2(  64), ceillog2(8)> dig_rev_config_R8_N64_SW8   = { DIG_REV_PERM_CONFIG_N64_SW8_R8       };
        perm_config<DIGIT_REV_NUM_STAGE_N512_SW8_R8 , 8, ceillog2( 512), ceillog2(8)> dig_rev_config_R8_N512_SW8  = { DIG_REV_PERM_CONFIG_N512_SW8_R8      };
        perm_config<DIGIT_REV_NUM_STAGE_N4096_SW8_R8, 8, ceillog2(4096), ceillog2(8)> dig_rev_config_R8_N4096_SW8 = { DIG_REV_PERM_CONFIG_N4096_SW8_R8     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW2_R2  , 2, ceillog2(  16), ceillog2(2)> stride_config_R2_N16_SW2   = { STRIDE_PERM_CONFIG_N16_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW2_R2  , 2, ceillog2(  32), ceillog2(2)> stride_config_R2_N32_SW2   = { STRIDE_PERM_CONFIG_N32_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW2_R2  , 2, ceillog2(  64), ceillog2(2)> stride_config_R2_N64_SW2   = { STRIDE_PERM_CONFIG_N64_SW2_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW2_R2 , 2, ceillog2( 128), ceillog2(2)> stride_config_R2_N128_SW2  = { STRIDE_PERM_CONFIG_N128_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW2_R2 , 2, ceillog2( 256), ceillog2(2)> stride_config_R2_N256_SW2  = { STRIDE_PERM_CONFIG_N256_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW2_R2 , 2, ceillog2( 512), ceillog2(2)> stride_config_R2_N512_SW2  = { STRIDE_PERM_CONFIG_N512_SW2_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW2_R2, 2, ceillog2(1024), ceillog2(2)> stride_config_R2_N1024_SW2 = { STRIDE_PERM_CONFIG_N1024_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW2_R2, 2, ceillog2(2048), ceillog2(2)> stride_config_R2_N2048_SW2 = { STRIDE_PERM_CONFIG_N2048_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW2_R2, 2, ceillog2(4096), ceillog2(2)> stride_config_R2_N4096_SW2 = { STRIDE_PERM_CONFIG_N4096_SW2_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R2  , 4, ceillog2(  16), ceillog2(4)> stride_config_R2_N16_SW4   = { STRIDE_PERM_CONFIG_N16_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW4_R2  , 4, ceillog2(  32), ceillog2(4)> stride_config_R2_N32_SW4   = { STRIDE_PERM_CONFIG_N32_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R2  , 4, ceillog2(  64), ceillog2(4)> stride_config_R2_N64_SW4   = { STRIDE_PERM_CONFIG_N64_SW4_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW4_R2 , 4, ceillog2( 128), ceillog2(4)> stride_config_R2_N128_SW4  = { STRIDE_PERM_CONFIG_N128_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R2 , 4, ceillog2( 256), ceillog2(4)> stride_config_R2_N256_SW4  = { STRIDE_PERM_CONFIG_N256_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW4_R2 , 4, ceillog2( 512), ceillog2(4)> stride_config_R2_N512_SW4  = { STRIDE_PERM_CONFIG_N512_SW4_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R2, 4, ceillog2(1024), ceillog2(4)> stride_config_R2_N1024_SW4 = { STRIDE_PERM_CONFIG_N1024_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW4_R2, 4, ceillog2(2048), ceillog2(4)> stride_config_R2_N2048_SW4 = { STRIDE_PERM_CONFIG_N2048_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R2, 4, ceillog2(4096), ceillog2(4)> stride_config_R2_N4096_SW4 = { STRIDE_PERM_CONFIG_N4096_SW4_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R2  , 8, ceillog2(  16), ceillog2(8)> stride_config_R2_N16_SW8   = { STRIDE_PERM_CONFIG_N16_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW8_R2  , 8, ceillog2(  32), ceillog2(8)> stride_config_R2_N32_SW8   = { STRIDE_PERM_CONFIG_N32_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R2  , 8, ceillog2(  64), ceillog2(8)> stride_config_R2_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R2       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW8_R2 , 8, ceillog2( 128), ceillog2(8)> stride_config_R2_N128_SW8  = { STRIDE_PERM_CONFIG_N128_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R2 , 8, ceillog2( 256), ceillog2(8)> stride_config_R2_N256_SW8  = { STRIDE_PERM_CONFIG_N256_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R2 , 8, ceillog2( 512), ceillog2(8)> stride_config_R2_N512_SW8  = { STRIDE_PERM_CONFIG_N512_SW8_R2      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R2, 8, ceillog2(1024), ceillog2(8)> stride_config_R2_N1024_SW8 = { STRIDE_PERM_CONFIG_N1024_SW8_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW8_R2, 8, ceillog2(2048), ceillog2(8)> stride_config_R2_N2048_SW8 = { STRIDE_PERM_CONFIG_N2048_SW8_R2     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R2, 8, ceillog2(4096), ceillog2(8)> stride_config_R2_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R2     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R4  , 4, ceillog2(  16), ceillog2(4)> stride_config_R4_N16_SW4   = { STRIDE_PERM_CONFIG_N16_SW4_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R4  , 4, ceillog2(  64), ceillog2(4)> stride_config_R4_N64_SW4   = { STRIDE_PERM_CONFIG_N64_SW4_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R4 , 4, ceillog2( 256), ceillog2(4)> stride_config_R4_N256_SW4  = { STRIDE_PERM_CONFIG_N256_SW4_R4      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R4, 4, ceillog2(1024), ceillog2(4)> stride_config_R4_N1024_SW4 = { STRIDE_PERM_CONFIG_N1024_SW4_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R4, 4, ceillog2(4096), ceillog2(4)> stride_config_R4_N4096_SW4 = { STRIDE_PERM_CONFIG_N4096_SW4_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R4  , 8, ceillog2(  16), ceillog2(8)> stride_config_R4_N16_SW8   = { STRIDE_PERM_CONFIG_N16_SW8_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R4  , 8, ceillog2(  64), ceillog2(8)> stride_config_R4_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R4       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R4 , 8, ceillog2( 256), ceillog2(8)> stride_config_R4_N256_SW8  = { STRIDE_PERM_CONFIG_N256_SW8_R4      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R4, 8, ceillog2(1024), ceillog2(8)> stride_config_R4_N1024_SW8 = { STRIDE_PERM_CONFIG_N1024_SW8_R4     };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R4, 8, ceillog2(4096), ceillog2(8)> stride_config_R4_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R4     };

        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R8  , 8, ceillog2(  64), ceillog2(8)> stride_config_R8_N64_SW8   = { STRIDE_PERM_CONFIG_N64_SW8_R8       };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R8 , 8, ceillog2( 512), ceillog2(8)> stride_config_R8_N512_SW8  = { STRIDE_PERM_CONFIG_N512_SW8_R8      };
        perm_config<STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R8, 8, ceillog2(4096), ceillog2(8)> stride_config_R8_N4096_SW8 = { STRIDE_PERM_CONFIG_N4096_SW8_R8     };
        
        if(RADIX == 2){
            if(SW == 2){
                switch (SIZE) {
                    case   16: sbg_radix_fft<dType, tType, iType,   16, 2, 2, GS, ceillog2(  16), ceillog2(2), DIGIT_REV_NUM_STAGE_N16_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW2_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N16_SW2  , stride_config_R2_N16_SW2  , Q_twiddle-3); break;
                    case   32: sbg_radix_fft<dType, tType, iType,   32, 2, 2, GS, ceillog2(  32), ceillog2(2), DIGIT_REV_NUM_STAGE_N32_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW2_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N32_SW2  , stride_config_R2_N32_SW2  , Q_twiddle-3); break;
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 2, 2, GS, ceillog2(  64), ceillog2(2), DIGIT_REV_NUM_STAGE_N64_SW2_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW2_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N64_SW2  , stride_config_R2_N64_SW2  , Q_twiddle-3); break;
                    case  128: sbg_radix_fft<dType, tType, iType,  128, 2, 2, GS, ceillog2( 128), ceillog2(2), DIGIT_REV_NUM_STAGE_N128_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW2_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N128_SW2 , stride_config_R2_N128_SW2 , Q_twiddle-3); break;
                    case  256: sbg_radix_fft<dType, tType, iType,  256, 2, 2, GS, ceillog2( 256), ceillog2(2), DIGIT_REV_NUM_STAGE_N256_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW2_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N256_SW2 , stride_config_R2_N256_SW2 , Q_twiddle-3); break;
                    case  512: sbg_radix_fft<dType, tType, iType,  512, 2, 2, GS, ceillog2( 512), ceillog2(2), DIGIT_REV_NUM_STAGE_N512_SW2_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW2_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N512_SW2 , stride_config_R2_N512_SW2 , Q_twiddle-3); break;
                    case 1024: sbg_radix_fft<dType, tType, iType, 1024, 2, 2, GS, ceillog2(1024), ceillog2(2), DIGIT_REV_NUM_STAGE_N1024_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW2_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N1024_SW2, stride_config_R2_N1024_SW2, Q_twiddle-3); break;
                    case 2048: sbg_radix_fft<dType, tType, iType, 2048, 2, 2, GS, ceillog2(2048), ceillog2(2), DIGIT_REV_NUM_STAGE_N2048_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW2_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N2048_SW2, stride_config_R2_N2048_SW2, Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 2, 2, GS, ceillog2(4096), ceillog2(2), DIGIT_REV_NUM_STAGE_N4096_SW2_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW2_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N4096_SW2, stride_config_R2_N4096_SW2, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }else if(SW == 4){
                switch (SIZE) {
                    case   16: sbg_radix_fft<dType, tType, iType,   16, 2, 4, GS, ceillog2(  16), ceillog2(4), DIGIT_REV_NUM_STAGE_N16_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N16_SW4  , stride_config_R2_N16_SW4  , Q_twiddle-3); break;
                    case   32: sbg_radix_fft<dType, tType, iType,   32, 2, 4, GS, ceillog2(  32), ceillog2(4), DIGIT_REV_NUM_STAGE_N32_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW4_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N32_SW4  , stride_config_R2_N32_SW4  , Q_twiddle-3); break;
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 2, 4, GS, ceillog2(  64), ceillog2(4), DIGIT_REV_NUM_STAGE_N64_SW4_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N64_SW4  , stride_config_R2_N64_SW4  , Q_twiddle-3); break;
                    case  128: sbg_radix_fft<dType, tType, iType,  128, 2, 4, GS, ceillog2( 128), ceillog2(4), DIGIT_REV_NUM_STAGE_N128_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW4_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N128_SW4 , stride_config_R2_N128_SW4 , Q_twiddle-3); break;
                    case  256: sbg_radix_fft<dType, tType, iType,  256, 2, 4, GS, ceillog2( 256), ceillog2(4), DIGIT_REV_NUM_STAGE_N256_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N256_SW4 , stride_config_R2_N256_SW4 , Q_twiddle-3); break;
                    case  512: sbg_radix_fft<dType, tType, iType,  512, 2, 4, GS, ceillog2( 512), ceillog2(4), DIGIT_REV_NUM_STAGE_N512_SW4_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW4_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N512_SW4 , stride_config_R2_N512_SW4 , Q_twiddle-3); break;
                    case 1024: sbg_radix_fft<dType, tType, iType, 1024, 2, 4, GS, ceillog2(1024), ceillog2(4), DIGIT_REV_NUM_STAGE_N1024_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N1024_SW4, stride_config_R2_N1024_SW4, Q_twiddle-3); break;
                    case 2048: sbg_radix_fft<dType, tType, iType, 2048, 2, 4, GS, ceillog2(2048), ceillog2(4), DIGIT_REV_NUM_STAGE_N2048_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW4_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N2048_SW4, stride_config_R2_N2048_SW4, Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 2, 4, GS, ceillog2(4096), ceillog2(4), DIGIT_REV_NUM_STAGE_N4096_SW4_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N4096_SW4, stride_config_R2_N4096_SW4, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            } else if(SW == 8){
                switch (SIZE) {
                    case   16: sbg_radix_fft<dType, tType, iType,   16, 2, 8, GS, ceillog2(  16), ceillog2(8), DIGIT_REV_NUM_STAGE_N16_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N16_SW8  , stride_config_R2_N16_SW8  , Q_twiddle-3); break;
                    case   32: sbg_radix_fft<dType, tType, iType,   32, 2, 8, GS, ceillog2(  32), ceillog2(8), DIGIT_REV_NUM_STAGE_N32_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N32_SW8_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N32_SW8  , stride_config_R2_N32_SW8  , Q_twiddle-3); break;
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 2, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R2  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R2  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N64_SW8  , stride_config_R2_N64_SW8  , Q_twiddle-3); break;
                    case  128: sbg_radix_fft<dType, tType, iType,  128, 2, 8, GS, ceillog2( 128), ceillog2(8), DIGIT_REV_NUM_STAGE_N128_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N128_SW8_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N128_SW8 , stride_config_R2_N128_SW8 , Q_twiddle-3); break;
                    case  256: sbg_radix_fft<dType, tType, iType,  256, 2, 8, GS, ceillog2( 256), ceillog2(8), DIGIT_REV_NUM_STAGE_N256_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N256_SW8 , stride_config_R2_N256_SW8 , Q_twiddle-3); break;
                    case  512: sbg_radix_fft<dType, tType, iType,  512, 2, 8, GS, ceillog2( 512), ceillog2(8), DIGIT_REV_NUM_STAGE_N512_SW8_R2 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R2 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N512_SW8 , stride_config_R2_N512_SW8 , Q_twiddle-3); break;
                    case 1024: sbg_radix_fft<dType, tType, iType, 1024, 2, 8, GS, ceillog2(1024), ceillog2(8), DIGIT_REV_NUM_STAGE_N1024_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N1024_SW8, stride_config_R2_N1024_SW8, Q_twiddle-3); break;
                    case 2048: sbg_radix_fft<dType, tType, iType, 2048, 2, 8, GS, ceillog2(2048), ceillog2(8), DIGIT_REV_NUM_STAGE_N2048_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N2048_SW8_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N2048_SW8, stride_config_R2_N2048_SW8, Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 2, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R2, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R2>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R2_N4096_SW8, stride_config_R2_N4096_SW8, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }else if(RADIX == 4){
            if(SW == 4){
                switch (SIZE) {
                    case   16: sbg_radix_fft<dType, tType, iType,   16, 4, 4, GS, ceillog2(  16), ceillog2(4), DIGIT_REV_NUM_STAGE_N16_SW4_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW4_R4  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N16_SW4  , stride_config_R4_N16_SW4  , Q_twiddle-3); break;
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 4, 4, GS, ceillog2(  64), ceillog2(4), DIGIT_REV_NUM_STAGE_N64_SW4_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW4_R4  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N64_SW4  , stride_config_R4_N64_SW4  , Q_twiddle-3); break;
                    case  256: sbg_radix_fft<dType, tType, iType,  256, 4, 4, GS, ceillog2( 256), ceillog2(4), DIGIT_REV_NUM_STAGE_N256_SW4_R4 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW4_R4 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N256_SW4 , stride_config_R4_N256_SW4 , Q_twiddle-3); break;
                    case 1024: sbg_radix_fft<dType, tType, iType, 1024, 4, 4, GS, ceillog2(1024), ceillog2(4), DIGIT_REV_NUM_STAGE_N1024_SW4_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW4_R4>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N1024_SW4, stride_config_R4_N1024_SW4, Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 4, 4, GS, ceillog2(4096), ceillog2(4), DIGIT_REV_NUM_STAGE_N4096_SW4_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW4_R4>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N4096_SW4, stride_config_R4_N4096_SW4, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            } else if(SW == 8){
                switch (SIZE) {
                    case   16: sbg_radix_fft<dType, tType, iType,   16, 4, 8, GS, ceillog2(  16), ceillog2(8), DIGIT_REV_NUM_STAGE_N16_SW8_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N16_SW8_R4  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N16_SW8  , stride_config_R4_N16_SW8  , Q_twiddle-3); break;
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 4, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R4  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R4  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N64_SW8  , stride_config_R4_N64_SW8  , Q_twiddle-3); break;
                    case  256: sbg_radix_fft<dType, tType, iType,  256, 4, 8, GS, ceillog2( 256), ceillog2(8), DIGIT_REV_NUM_STAGE_N256_SW8_R4 , STRIDE_PERM_SWITCH_NUM_STAGE_N256_SW8_R4 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N256_SW8 , stride_config_R4_N256_SW8 , Q_twiddle-3); break;
                    case 1024: sbg_radix_fft<dType, tType, iType, 1024, 4, 8, GS, ceillog2(1024), ceillog2(8), DIGIT_REV_NUM_STAGE_N1024_SW8_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N1024_SW8_R4>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N1024_SW8, stride_config_R4_N1024_SW8, Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 4, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R4, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R4>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R4_N4096_SW8, stride_config_R4_N4096_SW8, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }else if(RADIX == 8){
            if(SW == 8){
                switch (SIZE) {
                    case   64: sbg_radix_fft<dType, tType, iType,   64, 8, 8, GS, ceillog2(  64), ceillog2(8), DIGIT_REV_NUM_STAGE_N64_SW8_R8  , STRIDE_PERM_SWITCH_NUM_STAGE_N64_SW8_R8  >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R8_N64_SW8  , stride_config_R8_N64_SW8  , Q_twiddle-3); break;
                    case  512: sbg_radix_fft<dType, tType, iType,  512, 8, 8, GS, ceillog2( 512), ceillog2(8), DIGIT_REV_NUM_STAGE_N512_SW8_R8 , STRIDE_PERM_SWITCH_NUM_STAGE_N512_SW8_R8 >(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R8_N512_SW8 , stride_config_R8_N512_SW8 , Q_twiddle-3); break;
                    case 4096: sbg_radix_fft<dType, tType, iType, 4096, 8, 8, GS, ceillog2(4096), ceillog2(8), DIGIT_REV_NUM_STAGE_N4096_SW8_R8, STRIDE_PERM_SWITCH_NUM_STAGE_N4096_SW8_R8>(dinI, dinQ, outI, outQ, tc, ts, dig_rev_config_R8_N4096_SW8, stride_config_R8_N4096_SW8, Q_twiddle-3); break;
                    default:
                        exit( EXIT_FAILURE );
                }
            }
        }

        for(int i = 0; i < SIZE; i++){
            data.I.data()[i] = ((double) outI[i]) / (((dType) 1) << Q_in);
            data.Q.data()[i] = ((double) outQ[i]) / (((dType) 1) << Q_in);
        }

        delete[] dinI;
        delete[] dinQ;
        delete[] outI;
        delete[] outQ;
    }
    else{
        std::cout << "(EE) Error in lib_fft_perm function, the FFT implementation does not exist !" << std::endl;
        std::cout << "(EE) q_input    =  [" << q_input    << "]" << std::endl;
        std::cout << "(EE) q_internal =  [" << q_internal << "]" << std::endl;
        std::cout << "(EE) q_rom      =  [" << q_rom      << "]" << std::endl;
        std::cout << "(EE) This configuration is not supported by the FFT implementation model =  [" << p.toString("model") << "]." << std::endl;
        exit( EXIT_FAILURE );
    }
}