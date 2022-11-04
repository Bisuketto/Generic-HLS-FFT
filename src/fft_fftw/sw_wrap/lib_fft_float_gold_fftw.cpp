/*
*****************************************
* @author       ALMORIN Hugues
*
* @date         2021-05-10
* @brief        FFT C Module
* @details      This is a C module for FFT
*
* @version      Revision 0.0.0 - File Created
*
* @copyright    Copyright ARELIS 2020
*****************************************
*/

/***************************** Include Files *********************************/

#include "lib_fft_float_gold_fftw.hpp"

/************************** Constant Definitions *****************************/

/**************************** Type Definitions *******************************/

/***************** Macros (Inline Functions) Definitions *********************/

/************************** Function Prototypes ******************************/

/************************** Variable Definitions *****************************/

/*****************************************************************************/

void lib_fft_float_gold_fftw(DataVector& data, Parameters& p){

    const int32_t FFT_SIZE = p.toInt("FFT_SIZE");
    const std::string q_input    = p.toString("q_input");

    std::vector<double> samples_i(FFT_SIZE);
    std::vector<double> samples_q(FFT_SIZE);
    std::vector<double> fft_i(FFT_SIZE);
    std::vector<double> fft_q(FFT_SIZE);

    if((q_input == "double") || (q_input == "float")){
        for(int i = 0; i < FFT_SIZE; i++){
            samples_i[i] = data.I[i];
            samples_q[i] = data.Q[i];//samples_in_Q.read();
        }
    }
    else{
        const uint8_t Q_in = p.toInt("q_input");
        for(int i = 0; i < FFT_SIZE; i++){
            samples_i[i] = std::round(data.I[i]*(1 << (Q_in-1)));
            samples_q[i] = std::round(data.Q[i]*(1 << (Q_in-1)));//samples_in_Q.read();
        }
    }

    fft_fftw3 fft(FFT_SIZE);
    fft.load(0, samples_i, samples_q);
    fft.process();
    fft.store(fft_i, fft_q);
    
    if((q_input == "double") || (q_input == "float")){
        for(int i = 0; i < FFT_SIZE; i++){
            data.I[i]=fft_i[i];
            data.Q[i]=fft_q[i];
        }
    }
    else{    
        const uint8_t Q_in = p.toInt("q_input");
        for(int i = 0; i < FFT_SIZE; i++){
            data.I[i]=fft_i[i]/(1<<(Q_in-1));
            data.Q[i]=fft_q[i]/(1<<(Q_in-1));
        }
    }
    
}