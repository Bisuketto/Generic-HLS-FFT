/**
 * @file SIG_gen.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This tool generates examples signals for testing the FFT
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
#include <cstdint>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char* argv[])
{
    printf("# FFT input signal generatordevelopped by Bertrand LE GAL and Hugues ALMORIN (2022)\n");
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

    if(argc < 2){
        std::cout << "Usage : SIG_gen [dest_dir]" << std::endl;
        return -1;
    }
    std::string outdir(argv[1]);

    std::cout << "# Writing signals to " << outdir << std::endl;

    printf("#\n");

    const double pi = 3.1415926535897932384626433832795;

    const int32_t FFT_SIZE_MIN  =   8;
    const int32_t FFT_SIZE_MAX  = 65536;

    const int32_t NUM_BITS  = 16;
    const int32_t FRAC_BITS = NUM_BITS - 2;

    int result = mkdir(outdir.c_str(), 0777);

    printf("# Launching input signal generation\n");
    for(int N = FFT_SIZE_MIN; N <= FFT_SIZE_MAX; N *= 2)
    {
        printf("#  + FFT LENGTH IS : %4d\n", N);
        std::string n_direc  = outdir + "/N-" + std::to_string(N);

        // std::cout << "Creating " << n_direc;
        int result = mkdir(n_direc.c_str(), 0777);
        // if(result == 0)
        //     std::cout << " Done" << std::endl;
        // else{
        //     std::cout << " FAILED" << std::endl;
        //     return -1;
        // }

        std::string n_1sinus = n_direc + "/sinus-x1.txt";
        std::string n_2sinus = n_direc + "/sinus-x2.txt";
        std::string n_1ramp  = n_direc + "/ramp-x1.txt";
        std::string n_2ramp  = n_direc + "/ramp-x2.txt";
        std::string n_random = n_direc + "/noise.txt";

        FILE* f1 = fopen(n_1sinus.c_str(), "w");
        FILE* f2 = fopen(n_2sinus.c_str(),"w");
        FILE* f3 = fopen(n_1ramp.c_str(),"w");
        FILE* f4 = fopen(n_2ramp.c_str(),"w");
        FILE* f5 = fopen(n_random.c_str(),"w");

        for(int i = 0; i < N; i += 1)
        {
            const double anglex1 = (2*pi*i)/(double)N;
            const double anglex2 = (4*pi*i)/(double)N;
            const double sinx1   = sin(anglex1);
            const double sinx2   = sin(anglex2);
            const double rampx1  = (double) i      /(double) (N-1);
            const double rampx2  = (double)(i%(N/2)/(double)((N/2)-1));
            const double noise   = (double)rand()/(double)RAND_MAX;

            fprintf(f1, "%d %+23.23f %+23.23f\n", i, sinx1,  0.0f);
            fprintf(f2, "%d %+23.23f %+23.23f\n", i, sinx2,  0.0f);
            fprintf(f3, "%d %+23.23f %+23.23f\n", i, rampx1, 0.0f);
            fprintf(f4, "%d %+23.23f %+23.23f\n", i, rampx2, 0.0f);
            fprintf(f5, "%d %+23.23f %+23.23f\n", i, noise,  0.0f);
        }
        fclose( f1 );
        fclose( f2 );
        fclose( f3 );
        fclose( f4 );
        fclose( f5 );
    }
    printf("# Ending input signal generation\n");

    return EXIT_SUCCESS;
}
