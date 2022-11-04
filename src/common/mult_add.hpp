/**
 * @file mult_add.hpp
 * @author Hugues Almorin (hugues.almorin@arelis.com)
 * @brief This file contains mutiplier-additioner utility functions
 * @version 0.0.0
 * @date 2021-06-15
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */


#ifndef MULT_ADD_HPP_
#define MULT_ADD_HPP_

#include <cstdint>

/**
 * @brief Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam TW_TYPE Twiddle factor type
 * @tparam INT_TYPE Output type
 * @param a Input data
 * @param b Input tw
 * @param Q binary right shift for fixed point operation (ignored for floating point)
 * @return INT_TYPE a*b
 */
template<class DTYPE, class TW_TYPE, class INT_TYPE>
INT_TYPE mult(DTYPE  a, TW_TYPE b, uint8_t Q)
{
    return (a * b) >> Q;
}

/**
 * @brief Floating point Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam TW_TYPE Twiddle factor type
 * @tparam INT_TYPE Output type
 * @param a Input data
 * @param b Input tw
 * @param Q binary right shift for fixed point operation (ignored for floating point)
 * @return INT_TYPE a*b
 */
template< >
float mult(float a, float b, uint8_t Q);

/**
 * @brief Floating point Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam TW_TYPE Twiddle factor type
 * @tparam INT_TYPE Output type
 * @param a Input data
 * @param b Input tw
 * @param Q binary right shift for fixed point operation (ignored for floating point)
 * @return INT_TYPE a*b
 */
template< >
double mult(double a, double b, uint8_t Q);

/**
 * @brief Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Ouput
 * @param tw floating point litteral
 * @param b Input
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE tw*b
 */
template<class DTYPE, class INT_TYPE>
INT_TYPE mult_litt_f(float  tw, DTYPE b, uint8_t Q)
{
    const DTYPE fp_tw = tw * (1 << Q);
    return ((fp_tw * b) >> Q);
}

/**
 * @brief Floating point Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Ouput
 * @param tw floating point litteral
 * @param b Input
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE tw*b
 */
template< >
float mult_litt_f(float tw, float b, uint8_t Q);

/**
 * @brief Floating point Multiplier
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Ouput
 * @param tw floating point litteral
 * @param b Input
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE tw*b
 */
template< >
double mult_litt_f(float tw, double b, uint8_t Q);

/**
 * @brief Multiplier-additionner
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b+c*d
 */
template<class DTYPE, class TW_TYPE, class INT_TYPE>
INT_TYPE mult_add(DTYPE  a, TW_TYPE b, DTYPE c, TW_TYPE d, uint8_t Q)
{
    int k = 1;
    return (mult<DTYPE, TW_TYPE, INT_TYPE>(a,b,Q-k)+mult<DTYPE, TW_TYPE, INT_TYPE>(c,d,Q-k))>>k;                                                                   
}

/**
 * @brief Floating point Multiplier-additionner
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b+c*d
 */
template< >
float mult_add(float a, float b, float c, float d, uint8_t Q);

/**
 * @brief Floating point Multiplier-additionner
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b+c*d
 */
template< >
double mult_add(double a, double b, double c, double d, uint8_t Q);

/**
 * @brief Multiplier-substractor
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b-c*d
 */
template<class DTYPE, class TW_TYPE, class INT_TYPE>
INT_TYPE mult_sub(DTYPE  a, TW_TYPE b, DTYPE c, TW_TYPE d, uint8_t Q)
{
    int k = 1;
    return (mult<DTYPE, TW_TYPE, INT_TYPE>(a,b,Q-k)-mult<DTYPE, TW_TYPE, INT_TYPE>(c,d,Q-k))>>k;
}

/**
 * @brief Floating point Multiplier-substractor
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b-c*d
 */
template< >
float mult_sub(float a, float b, float c, float d, uint8_t Q);

/**
 * @brief Floating point Multiplier-substractor
 * 
 * @tparam DTYPE data type
 * @tparam TW_TYPE twiddle factor type
 * @tparam INT_TYPE output type
 * @param a data in a
 * @param b tw in b
 * @param c data in c
 * @param d tw in d
 * @param Q binary shift for fixed point operation
 * @return INT_TYPE a*b-c*d
 */
template< >
double mult_sub(double a, double b, double c, double d, uint8_t Q);

/**
 * @brief Scaling utility function
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Output type
 * @param in Input data
 * @return DTYPE scaled value
 */
template<class DTYPE, class INT_TYPE>
DTYPE util_scale(INT_TYPE in){
    return (in >> 1);
}

/**
 * @brief Floating point Scaling utility function
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Output type
 * @param in Input data
 * @return DTYPE scaled value
 */
template< >
float util_scale(float in);

/**
 * @brief Floating point Scaling utility function
 * 
 * @tparam DTYPE Data type
 * @tparam INT_TYPE Output type
 * @param in Input data
 * @return DTYPE scaled value
 */
template< >
double util_scale(double in);

/**
 * @brief Commplex scaling utility function
 * 
 * @tparam itype input type
 * @tparam stype scaling factor type
 * @tparam otype output type
 * @param in_cplx_I real input
 * @param in_cplx_Q imag input
 * @param out_cplx_I real output
 * @param out_cplx_Q imag output 
 * @param factor scaling factor
 */
template<class itype, class stype, class otype>
void scale_cplx(itype in_cplx_I, itype in_cplx_Q, otype& out_cplx_I, otype& out_cplx_Q, stype factor){
    out_cplx_I = in_cplx_I * factor;
    out_cplx_Q = in_cplx_Q * factor;
}

/**
 * @brief Complex multiplier, using 3 multiplier and 5 additioners.
 * 
 * @tparam atype First complex value input type
 * @tparam btype Second complex value input type
 * @tparam otype Complex value output type
 * @param a_r First complex value input data (real)
 * @param a_i First complex value input data (imag)
 * @param b_r Second complex value input data (real)
 * @param b_i Second complex value input data (imag)
 * @param o_r Complex value output data (real)
 * @param o_i Complex value output data (imag)
 * @param Q quantization factor
 */
template<class atype, class btype, class otype, class itype>
void cmult(atype a_r, atype a_i, btype b_r, btype b_i, otype& o_r, otype& o_i, int Q){
    itype x = (a_r + a_i) >> 1;
    itype y = (b_i - b_r) >> 1;
    itype z = (b_r + b_i) >> 1;
    itype k1 = (b_r * x) >> (Q-1);
    itype k2 = (a_r * y) >> (Q-1);
    itype k3 = (a_i * z) >> (Q-1);
    o_r = (k1 - k3);
    o_i = (k1 + k2);
}

/**
 * @brief Floating point Complex multiplier, using 3 multiplier and 5 additioners.
 * 
 * @tparam atype First complex value input type
 * @tparam btype Second complex value input type
 * @tparam otype Complex value output type
 * @param a_r First complex value input data (real)
 * @param a_i First complex value input data (imag)
 * @param b_r Second complex value input data (real)
 * @param b_i Second complex value input data (imag)
 * @param o_r Complex value output data (real)
 * @param o_i Complex value output data (imag)
 * @param Q quantization factor
 */
template< >
void cmult<float, float, float, float>(float a_r, float a_i, float b_r, float b_i, float& o_r, float& o_i, int Q);

/**
 * @brief Floating point Complex multiplier, using 3 multiplier and 5 additioners.
 * 
 * @tparam atype First complex value input type
 * @tparam btype Second complex value input type
 * @tparam otype Complex value output type
 * @param a_r First complex value input data (real)
 * @param a_i First complex value input data (imag)
 * @param b_r Second complex value input data (real)
 * @param b_i Second complex value input data (imag)
 * @param o_r Complex value output data (real)
 * @param o_i Complex value output data (imag)
 * @param Q quantization factor
 */
template< >
void cmult<double, double, double, double>(double a_r, double a_i, double b_r, double b_i, double& o_r, double& o_i, int Q);

#endif // MULT_ADD_HPP_