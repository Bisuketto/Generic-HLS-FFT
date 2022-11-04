/**
 * @file fft_fftw3.hpp
 * @author Lamia DACHI, Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains functions using the fftw3 library
 * @version 0.0.0
 * @date 2020-07-06
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef TX_CHAIN_FFT_FFTW3_H
#define TX_CHAIN_FFT_FFTW3_H

#include <cstdint>
#include <vector>
#include <fftw3.h>

using namespace std;

class fft_fftw3 {
public :
    fft_fftw3(const uint32_t _M);
    ~fft_fftw3();

    void load(uint32_t ech_exam, vector<double>& _inR, vector<double>& _inI);
    void process();
    void store(vector<double>& _outR, vector<double>& _outI);

private :
    uint32_t  length;
    fftw_plan p;

    fftw_complex* in;
    fftw_complex* out;
};

#endif //TX_CHAIN_FFT_FFTW3_H
