/**
 * @file LUT_gen.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This tools generates LUT files for the twiddle factors
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


unsigned int reverse_bits(unsigned int input, int M) {
	int i, rev = 0;
	rev_bits_loop : for (i = 0; i < M; i++) {
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

void index_permut(
    unsigned int R,
    unsigned int N,
    unsigned int* index)
{
    unsigned int m = R;
    unsigned int n = N/R;

    for(unsigned int i = 0; i < m; i++){
        for(unsigned j = 0; j < n; j++){
            index[m*j+i] = n*i+j;
        }
    }
}

int main(int argc, char* argv[])
{
    printf("# FFT ROM generation program developped by Bertrand LE GAL and Hugues ALMORIN (2022)\n");
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
        std::cout << "Usage : LUT_gen [dest_dir]" << std::endl;
        return -1;
    }
    std::string outdir(argv[1]);

    std::cout << "(II) - Writing roms to " << outdir << std::endl;


    printf("#\n");

    const double pi = 3.1415926535897932384626433832795;

    const int32_t FFT_SIZE_MIN  =     8;
    const int32_t FFT_SIZE_MAX  = 8192;

    const int32_t COEF_BIT_MIN  =  8;
    const int32_t COEF_BIT_MAX  = 32;


    int result = mkdir(outdir.c_str(), 0777);
    FILE* froms = fopen((outdir + "/roms.h").c_str(), "w");

    printf("# Launching ROM generation\n");
    printf("#\n");

    //
    // On genere toutes les tailles de FFTs demandée par l'utilisateur
    //

    fprintf(froms, "// COS LUT\n");
    printf("#   COSINUS LOOK-UP TABLES\n");
    for(int N = FFT_SIZE_MIN; N <= FFT_SIZE_MAX; N *= 2) {

        printf("#   + FFT LENGTH IS : %5d : ", N);
        std::string n_direc = outdir + "/N-" + std::to_string(N);
        std::string f_direc = n_direc + "/float";
        std::string cos_fl_lut = f_direc + "/cos_lut_N" + std::to_string(N) + "_float.h";
        fprintf(froms, "#include \"./N-%d/float/cos_lut_N%d_float.h\"\n", N, N);

        int result = mkdir(n_direc.c_str(), 0777);

        for (int NUM_BITS = COEF_BIT_MIN; NUM_BITS <= COEF_BIT_MAX; NUM_BITS += 1)
        {
            const int32_t FRAC_BITS = NUM_BITS - 2;
            printf("Q%d.%d ", 2, FRAC_BITS);

            std::string q_direc = n_direc + "/Q" + std::to_string(NUM_BITS);

            mkdir(q_direc.c_str(), 0777);
            mkdir(f_direc.c_str(), 0777);

            std::string cos_fixed  = q_direc + "/cos.txt";
            std::string cos_report = q_direc + "/cos_conv_report.txt";
            std::string cos_double = f_direc + "/cos.txt";
            std::string cos_fp_lut = q_direc + "/cos_lut_N" + std::to_string(N) + "_Q2" + "_" + std::to_string(FRAC_BITS) + ".h";
            fprintf(froms, "#include \"./N-%d/Q%d/cos_lut_N%d_Q2_%d.h\"\n", N, NUM_BITS, N, FRAC_BITS);

            FILE *fcos_fixed  = fopen(cos_fixed.c_str(),  "w");
            FILE *fcos_double = fopen(cos_double.c_str(), "w");
            FILE *fcos_fl_lut = fopen(cos_fl_lut.c_str(), "w");
            FILE *fcos_fp_lut = fopen(cos_fp_lut.c_str(), "w");
            FILE *fcos_report = fopen(cos_report.c_str(), "w");

            fprintf(fcos_fp_lut, "#define COS_LUT_N%d_Q%d_%d {", N, 2, FRAC_BITS);
            fprintf(fcos_fl_lut, "#define COS_LUT_N%d_float {", N);
            for (int i = 0; i < N; i += 1)
            {
                const double angle    = 2 * pi * i / (double) N;
                const double data_dp  = cos(angle);
                const int32_t data_fp = std::round(data_dp * ((double) (1 << (FRAC_BITS - 1))));
                const double data_sp  = ((double) data_fp) / ((double) (1 << (FRAC_BITS - 1)));
                const double data_er  = data_dp - data_sp;

                fprintf(fcos_fixed, "%d\n", data_fp);
                fprintf(fcos_fp_lut,"%d", data_fp);                
                fprintf(fcos_double, "%+23.23f\n", data_dp);
                fprintf(fcos_fl_lut, "%+23.23f", data_dp);
                fprintf(fcos_report,
                        "N = %d  |  double (%+23.23f)  |  double (%+23.23f)  |  fixed(%+10d)  |  value(%+23.23f) | error(%+23.23f)\n",
                        i, angle, data_dp, data_fp, data_sp, data_er);
                //            fprintf(f_report, "double (%+23.23f)  |  fixed(%+10d)  |  value(%+23.23f) | error(%+23.23f)\n", data_dp, data_fp, data_sp, data_er);
                //            fprintf("double (%+23.23f)  |  fixed(%+10d)  |  value(%+23.23f) | error(%+23.23f)\n", data_dp, data_fp, data_sp, data_er);
                if(i < N-1){
                    fprintf(fcos_fp_lut, ", ");
                    fprintf(fcos_fl_lut, ", ");
                }
            }
            fprintf(fcos_fp_lut, "}");
            fprintf(fcos_fl_lut, "}");
            fclose(fcos_fp_lut);
            fclose(fcos_fl_lut);
            fclose(fcos_fixed);
            fclose(fcos_double);
            fclose(fcos_report);
        }
        printf("\n");
    }
    printf("#\n");


    //
    // On genere toutes les tailles de FFTs demandée par l'utilisateur
    //
    fprintf(froms, "\n// SIN LUT\n");
    printf("#   SINUS LOOK-UP TABLES\n");
    for(int N = FFT_SIZE_MIN; N <= FFT_SIZE_MAX; N *= 2) {

        printf("#   + FFT LENGTH IS : %5d : ", N);
        std::string n_direc = outdir + "/N-" + std::to_string(N);
        std::string f_direc = n_direc + "/float";
        std::string sin_fl_lut = f_direc + "/sin_lut_N" + std::to_string(N) + "_float.h";
        fprintf(froms, "#include \"./N-%d/float/sin_lut_N%d_float.h\"\n", N, N);

        int result = mkdir(n_direc.c_str(), 0777);

        for (int NUM_BITS = COEF_BIT_MIN; NUM_BITS <= COEF_BIT_MAX; NUM_BITS += 1)
        {
            const int32_t FRAC_BITS = NUM_BITS - 2;
            printf("Q%d.%d ", 2, FRAC_BITS);

            std::string q_direc = n_direc + "/Q" + std::to_string(NUM_BITS);
            int result = mkdir(q_direc.c_str(), 0777);

            std::string sin_fixed  = q_direc + "/sin.txt";
            std::string sin_report = q_direc + "/sin_conv_report.txt";
            std::string sin_double = f_direc + "/sin.txt";
            std::string sin_fp_lut = q_direc + "/sin_lut_N" + std::to_string(N) + "_Q2" + "_" + std::to_string(FRAC_BITS) + ".h";
            fprintf(froms, "#include \"./N-%d/Q%d/sin_lut_N%d_Q2_%d.h\"\n", N, NUM_BITS, N, FRAC_BITS);

            FILE* f_fixed  = fopen(sin_fixed.c_str(), "w");
            FILE* f_double = fopen(sin_double.c_str(),"w");
            FILE* fsin_fl_lut = fopen(sin_fl_lut.c_str(),"w");
            FILE* fsin_fp_lut = fopen(sin_fp_lut.c_str(),"w");
            FILE* f_report = fopen(sin_report.c_str(),"w");

            fprintf(fsin_fp_lut, "#define SIN_LUT_N%d_Q%d_%d {", N, 2, FRAC_BITS);
            fprintf(fsin_fl_lut, "#define SIN_LUT_N%d_float {", N);
            for(int i = 0; i < N; i += 1)
            {
                const double  angle   = 2 * pi * i / (double)N;
                const double  data_dp = sin(angle);
                const int32_t data_fp = std::round(data_dp * ( (double) (1 << (FRAC_BITS-1)) ));
                const double  data_sp = ((double)data_fp)  / ( (double) (1 << (FRAC_BITS-1)) );
                const double  data_er = data_dp - data_sp;

                fprintf(f_fixed,  "%d\n", data_fp);
                fprintf(fsin_fp_lut,"%d", data_fp);                
                fprintf(f_double, "%+23.23f\n", data_dp);
                fprintf(fsin_fl_lut, "%+23.23f", data_dp);
                fprintf(f_report, "N = %d  |  double (%+23.23f)  |  double (%+23.23f)  |  fixed(%+10d)  |  value(%+23.23f) | error(%+23.23f)\n", i, angle, data_dp, data_fp, data_sp, data_er);
                if(i < N-1){
                    fprintf(fsin_fp_lut, ", ");
                    fprintf(fsin_fl_lut, ", ");
                }
            }
            fprintf(fsin_fp_lut, "}");
            fprintf(fsin_fl_lut, "}");
            fclose( fsin_fp_lut );
            fclose( fsin_fl_lut );
            fclose( f_fixed  );
            fclose( f_double );
            fclose( f_report );
        }
        printf("\n");
    }
    printf("#\n");

    fprintf(froms, "\n// BR LUT\n");
    printf("#   BIT-REVERSE LOOK-UP TABLES\n");
    for(int N = FFT_SIZE_MIN; N <= FFT_SIZE_MAX; N *= 2) {
        printf("#   + FFT LENGTH IS : %5d ", N);
        std::string n_direc = outdir + "/N-" + std::to_string(N);

        mkdir(n_direc.c_str(), 0777);
        
        std::string br_direc = n_direc + "/br";
        mkdir(br_direc.c_str(), 0777);

        std::string br_lut = br_direc + "/br_lut_N" + std::to_string(N) + ".h";
        fprintf(froms, "#include \"./N-%d/br/br_lut_N%d.h\"\n", N, N);

        FILE* fbr_lut = fopen(br_lut.c_str(), "w");

        fprintf(fbr_lut, "#define BR_LUT_N%d {", N);
        for (int i = 0; i < N; i++)
        {
            fprintf(fbr_lut, "%d", reverse_bits(i, log2(N)));

            if(i < N-1){
                fprintf(fbr_lut, ", ");
            }
        }
        fprintf(fbr_lut, "}");
        fclose(fbr_lut);

        printf("\n");
    }
    printf("#\n");

    fprintf(froms, "\n// PERM LUT\n");
    printf("#   RADIX-2 PERM STRIDE LOOK-UP TABLES\n");
    for(int N = FFT_SIZE_MIN; N <= FFT_SIZE_MAX; N *= 2) {
        printf("#   + FFT LENGTH IS : %5d ", N);
        std::string n_direc = outdir + "/N-" + std::to_string(N);

        mkdir(n_direc.c_str(), 0777);
        
        std::string perm_direc = n_direc + "/perm";
        mkdir(perm_direc.c_str(), 0777);

        std::string perm_lut = perm_direc + "/perm_lut_r2_N" + std::to_string(N) + ".h";
        fprintf(froms, "#include \"./N-%d/perm/perm_lut_r2_N%d.h\"\n", N, N);

        FILE* fperm_lut = fopen(perm_lut.c_str(), "w");

        unsigned int* permuted_index = new unsigned int[N];
        index_permut(2, N, permuted_index);

        fprintf(fperm_lut, "#define PERM_LUT_R2_N%d {", N);
        for (int i = 0; i < N; i++)
        {
            fprintf(fperm_lut, "%d", permuted_index[i]);

            if(i < N-1){
                fprintf(fperm_lut, ", ");
            }
        }
        fprintf(fperm_lut, "}");
        fclose(fperm_lut);
        delete[] permuted_index;

        printf("\n");
    }
    printf("#\n");
    
    fclose(froms);

    printf("# Ending ROM generation\n");

    return EXIT_SUCCESS;
}
