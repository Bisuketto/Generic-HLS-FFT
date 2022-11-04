/**
 * @file fft_sbg_radix_tb.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains the testbench for the FFT model.
 * @version 0.0.0
 * @date 2022-02-16
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include <iostream>
#include <cmath>
#include <complex>
#include "../../common/fft_fftw3.hpp"
#include "../module/fft_sbg_radix_hlstop.hpp"

template<class dType>
void fft_gold(dType* din_R, dType* din_I, dType* dout_R, dType* dout_I, int FFT_SIZE){

    std::vector<double> samples_i(FFT_SIZE);
    std::vector<double> samples_q(FFT_SIZE);
    std::vector<double> fft_i(FFT_SIZE);
    std::vector<double> fft_q(FFT_SIZE);

    
    for(int i = 0; i < FFT_SIZE; i++){
        samples_i[i] = din_R[i];
        samples_q[i] = din_I[i];
    }

    fft_fftw3 fft(FFT_SIZE);
    fft.load(0, samples_i, samples_q);
    fft.process();
    fft.store(fft_i, fft_q);
    
    for(int i = 0; i < FFT_SIZE; i++){
        dout_R[i]=fft_i[i];
        dout_I[i]=fft_q[i];
    }
}

int main(){
    std::cout << "FFT SBG RADIX TESTBENCH" << std::endl;
    SAMPLE_TYPE IN_R[NFFT], IN_I[NFFT];
    SAMPLE_TYPE OUT_R[NFFT], OUT_I[NFFT];

    double IN_gold_R[NFFT], IN_gold_I[NFFT];
    double OUT_gold_R[NFFT], OUT_gold_I[NFFT];

    for(int i = 0; i < NFFT; i++){
        IN_R[i] = sin(2*M_PI*i/NFFT)*(((SAMPLE_TYPE) 1) << (D_WIDTH-2));
        IN_I[i] = 0;
        IN_gold_R[i] = sin(2*M_PI*i/NFFT);
        IN_gold_I[i] = 0;
        // std::cout << IN_gold_R[i] << " " << IN_R[i] << std::endl;
    }

    fft_gold(IN_gold_R, IN_gold_I, OUT_gold_R, OUT_gold_I, NFFT);
    fft(IN_R, IN_I, OUT_R, OUT_I);

    double msle;
    std::complex<double> r;
    std::complex<double> m;
    
    for(int i = 0; i < NFFT; i++){
        r.real(OUT_gold_R[i]);
        r.imag(OUT_gold_I[i]);
        m.real(((double) OUT_R[i])/(((SAMPLE_TYPE) 1) << (D_WIDTH-2)));
        m.imag(((double) OUT_I[i])/(((SAMPLE_TYPE) 1) << (D_WIDTH-2)));
        msle += (1./NFFT)*pow(log10(std::abs(r)+1)-log10(std::abs(m)+1), 2);
        // printf("(%3d) GOLD : %+04.5lf %+04.5lf MEAS : %+04.5lf %+04.5lf\n", i, real(r), imag(r), real(m), imag(m));
        // std::cout << "N" << i << " GOLD R : " << real(r) << " I : " << imag(r) << " MEAS R : " << real(m) << " I : " << imag(m) << std::endl;
    }
    std::cout << "MSLE : " << msle << std::endl;

    return !(msle < 1e-5);
}
