/**
 * @file fft_calc.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This source contains a tool allowing to run the FFT model for selected parameters
 * @version 0.0.0
 * @date 2022-09-07
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include "../common/Parameters.hpp"
#include "../common/DataVector.hpp"

#include "../fft_fftw/sw_wrap/lib_fft_float_gold_fftw.hpp"
#include "../fft/sw_wrap/lib_fft_sbg_radix.hpp"


std::string replace(std::string str, std::string old, std::string news)
{
    int index = str.find(old);
    if (index == std::string::npos)
        return str;
    str.replace(index, old.size(), news);
    return str;
}

int main(int argc, char* argv[])
{
    printf("# FFT processing program developped by Bertrand LE GAL and Hugues ALMORIN (2022)\n");
    printf("#  + Binary generated : %s - %s\n", __DATE__, __TIME__);
#if defined(__clang__)
    /* Clang/LLVM. ---------------------------------------------- */
    printf("#  + Clang/LLVM version %d.%d.%d\n", __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    /* Intel ICC/ICPC. ------------------------------------------ */
    printf("#  + Intel ICC/ICPC version %d.%d\n", __INTEL_COMPILER, __INTEL_COMPILER_BUILD_DATE);
#elif defined(__GNUC__) || defined(__GNUG__)
	/* GNU GCC/G++. --------------------------------------------- */
    printf("#  + GNU GCC/G++ version %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER)
	/* Microsoft Visual Studio. --------------------------------- */
    printf("#  + Microsoft Visual Studio\n");
#else
    #error "#  + Undetected compiler !"
#endif
    printf("#\n");

    Parameters param;
    param.set("ifile",      "none");
    param.set("ofile",      "none");
    param.set("cfile",      "none");

    param.set("model",      "none");

    param.set("q_input",    "float");
    param.set("q_internal", "float");
    param.set("q_rom",      "float");
    param.set("q_rom_dir",  "float");
    param.set("SW",         "none");
    param.set("RADIX",      "none");


    for (uint32_t p = 1; p < argc; p++) {
        const std::string cmde = argv[p    ];
        const std::string arg1 = (p+1 < argc) ? argv[p + 1] : "";
        const std::string arg2 = (p+2 < argc) ? argv[p + 1] : "";
        const std::string arg3 = (p+3 < argc) ? argv[p + 1] : "";

        if (cmde == "--in-signal") {
            param.set("ifile", arg1);
            p += 1;
        } else if (cmde == "--out-signal") {
            param.set("ofile", arg1);
            p += 1;
        } else if (cmde == "--constants") {
            param.set("cfile", arg1);
            p += 1;
        } else if (cmde == "--input-quantif") {
            param.set("q_input", arg1);
            p += 1;
        } else if (cmde == "--internal-quantif") {
            param.set("q_internal", arg1);
            p += 1;
        } else if (cmde == "--rom-quantif") {
            param.set("q_rom", arg1);
            if( arg1 == "float" ) param.set("q_rom_dir", arg1);
            else                  param.set("q_rom_dir", "Q" + arg1);
            p += 1;
        } else if (cmde == "--fft-model") {
            param.set("model", arg1);
            p += 1;
        } else if (cmde == "--sw") {
            param.set("SW", arg1);
            p += 1;
        } else if (cmde == "--radix") {
            param.set("RADIX", arg1);
            p += 1;
        }else {
            printf("(EE) Unknown argument (%d) => [%s]\n", p, cmde.c_str());
            printf("(EE) Error in %s %d\n", __FILE__, __LINE__);
            exit(0);
        }
    }

    //
    // On charge les données d'entrée, cela permettra de connaitre la taille de
    // la fft. Cela permettra de générer le nom des fichiers ROM + output.
    //

    if(param.toString("ifile") == "none"){
        std::cout << "(EE) Please specify input file" << std::endl;
        return -1;
    }
    if(param.toString("ofile") == "none"){
        std::cout << "(EE) Please specify output file" << std::endl;
        return -1;
    }
    if(param.toString("cfile") == "none"){
        std::cout << "(EE) Please specify constants file" << std::endl;
        return -1;
    }
    if(param.toString("model") == "none"){
        std::cout << "(EE) Please specify FFT model" << std::endl;
        return -1;
    }

    if(param.toString("SW") == "none"){
        param.set("SW", "8");
    }
    if(param.toString("RADIX") == "none"){
        param.set("RADIX", "8");
    }

    DataVector  ii( param.toString("ifile") );

    const int32_t SIZE  = ii.I.size();
    const int32_t LOG2  = std::log2(SIZE);
    param.set("FFT_SIZE", SIZE);
    param.set("FFT_LOG2", LOG2);

    std::string cfile = param.toString("cfile");

    std::string ofile = param.toString("ofile");

    printf("#  + FFT configuration :\n");
    printf("#   - FFT model      : %s\n", param.toString("model").c_str());
    printf("#   - FFT size       : %s\n", param.toString("FFT_SIZE").c_str());
    printf("#   - FFT depth      : %s\n", param.toString("FFT_LOG2").c_str());
    printf("#   - I/O width      : %s\n", param.toString("q_input"   ).c_str());
    printf("#   - Internal width : %s\n", param.toString("q_internal").c_str());
    printf("#   - COS/SIN  width : %s\n", param.toString("q_rom"     ).c_str());
    printf("#\n");
    printf("#  + I/O file configuration :\n");
    printf("#   - LUT values   : %s\n", param.toString("cfile").c_str());
    printf("#   - Input signal : %s\n", param.toString("ifile").c_str());
    printf("#   - Ouput signal : %s\n", param.toString("ofile").c_str());
    printf("#\n");

    printf("#  > Launching FFT computation\n");

    const std::string fft_model = param.toString("model");
    if( fft_model == "sbg-radix" ) {
        lib_fft_sbg_radix_load_fft_roms( cfile, SIZE );
        lib_sbg_radix_fft(ii, param);
    }else if( fft_model == "fftw-gold" ) {
        lib_fft_float_gold_fftw(ii, param);
    }else{
        std::cout << "(EE) Error in main function, the FFT model does not exist !" << std::endl;
        std::cout << "(EE) model      =  [" << param.toString("model")      << "]" << std::endl;
        std::cout << "(EE) q_input    =  [" << param.toString("q_input")    << "]" << std::endl;
        std::cout << "(EE) q_internal =  [" << param.toString("q_internal") << "]" << std::endl;
        std::cout << "(EE) q_rom      =  [" << param.toString("q_rom")      << "]" << std::endl;
        exit( EXIT_FAILURE );
    }

    printf("#  > FFT computation ended\n");

    ii.save( ofile );

    printf("#  > FFT results saved\n");

    return 0;
}
