/**
 * @file mult_add.cpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief This file contains mutiplier-additioner utility functions
 * @version 0.0.0
 * @date 2021-06-29
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include "mult_add.hpp"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wunused-parameter"

template< >
float mult(float a, float b, uint8_t Q)
{
    return (a * b);    
}

template< >
double mult(double a, double b, uint8_t Q)
{
    return (a * b);    
}

template< >
float mult_litt_f(float tw, float b, uint8_t Q){
    return (tw * b);
}

template< >
double mult_litt_f(float tw, double b, uint8_t Q){
    return (tw * b);
}

template< >
float mult_add(float a, float b, float c, float d, uint8_t Q)
{
    return (a * b) + (c * d);    
}

template< >
double mult_add(double a, double b, double c, double d, uint8_t Q)
{
    return (a * b) + (c * d);    
}

template< >
float mult_sub(float a, float b, float c, float d, uint8_t Q)
{
    return (a * b) - (c * d);    
}

template< >
double mult_sub(double a, double b, double c, double d, uint8_t Q)
{
    return (a * b) - (c * d);
}

template< >
float util_scale(float in){
    return in;
}

template< >
double util_scale(double in){
    return in;
}

template< >
void cmult<float, float, float, float>(float a_r, float a_i, float b_r, float b_i, float& o_r, float& o_i, int Q){
    float x = a_r + a_i;
    float y = b_i - b_r;
    float z = b_r + b_i;
    float k1 = (b_r * x);
    float k2 = (a_r * y);
    float k3 = (a_i * z);
    o_r = (k1 - k3);
    o_i = (k1 + k2);
}

template< >
void cmult<double, double, double, double>(double a_r, double a_i, double b_r, double b_i, double& o_r, double& o_i, int Q){
    double x = a_r + a_i;
    double y = b_i - b_r;
    double z = b_r + b_i;
    double k1 = (b_r * x);
    double k2 = (a_r * y);
    double k3 = (a_i * z);
    o_r = (k1 - k3);
    o_i = (k1 + k2);
}

#pragma GCC diagnostic pop