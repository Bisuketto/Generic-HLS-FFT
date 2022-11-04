/**
 * @file fft_fftw3.cpp
 * @author Lamia DACHI, modified by Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains functions using the fftw3 library
 * @version 0.0.0
 * @date 2020-07-06
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */


#include "fft_fftw3.hpp"
fft_fftw3::fft_fftw3(const uint32_t _M)
{
    length = _M;
    in     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    out    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    p      = fftw_plan_dft_1d(length, in, out, FFTW_FORWARD, FFTW_MEASURE);  //FFTW_ESTIMATE
}

fft_fftw3::~fft_fftw3(){
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

void fft_fftw3::load(uint32_t ech_exam, vector<double>& _inR, vector<double>& _inI)
{
    for (uint32_t i=0; i < length; i++){
        in[i][0] = _inR[i + ech_exam];
        in[i][1] = _inI[i + ech_exam];
    }
}

void fft_fftw3::process()
{
    fftw_execute(p);
}

void fft_fftw3::store(vector<double>& _outR, vector<double>& _outI)
{
    for(uint32_t i=0; i < length; i++){
        _outR[i] = out[i][0];
        _outI[i] = out[i][1];
    }
}

