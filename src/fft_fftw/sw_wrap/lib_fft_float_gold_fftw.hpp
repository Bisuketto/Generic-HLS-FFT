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

#ifndef FFT_FFTW_H_
#define FFT_FFTW_H_

/***************************** Include Files *********************************/

#include "../../common/fft_fftw3.hpp"
#include "../../common/Parameters.hpp"
#include "../../common/DataVector.hpp"

/************************** Constant Definitions *****************************/

#define FFT_DEBUG //< Comment to hide debug printf's

#define FFT_VERSION_MAJOR 0
#define FFT_VERSION_MINOR 0
#define FFT_VERSION_PATCH 0

/**************************** Type Definitions *******************************/
 
// typedef std::complex<double> Complex;
// typedef std::valarray<Complex> CArray;


/***************** Macros (Inline Functions) Definitions *********************/

#ifdef FFT_DEBUG
    #define FFT_printf(...) printf(__VA_ARGS__)
#else
    #define FFT_printf(...)
#endif

/************************** Function Prototypes ******************************/

void lib_fft_float_gold_fftw(DataVector& data, Parameters& p);

/************************** Classes Definitions ******************************/

/*****************************************************************************/

#endif // FFT_FFTW_H_
