/**
 * @file const_logs.hpp
 * @author Hugues ALMORIN (hugues.almorin@arelis.com)
 * @brief log contexpr tools
 * @version 0.0.0
 * @date 2021-09-09
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef CONST_LOGS_HPP_
#define CONST_LOGS_HPP_

/**
 * @brief Compile time floorlog2
 * 
 * @param x Input
 * @return constexpr unsigned floor(log2(x))
 */
constexpr unsigned floorlog2(unsigned x)
{
    return x < 2 ? 0 : 1+floorlog2(x >> 1);
}

/**
 * @brief Compile time ceillog2
 * 
 * @param x Input
 * @return constexpr unsigned ceil(log2(x))
 */
constexpr unsigned ceillog2(unsigned x)
{
    return x < 2 ? 0 : floorlog2(x - 1) + 1;
}

/**
 * @brief Compile time floorlog4
 * 
 * @param x Input
 * @return constexpr unsigned floor(log4(x))
 */
constexpr unsigned floorlog4(unsigned x)
{
    return x < 4 ? 0 : 1+floorlog4(x >> 2);
}

/**
 * @brief Compile time ceillog4
 * 
 * @param x input
 * @return constexpr unsigned ceil(log4(x))
 */
constexpr unsigned ceillog4(unsigned x)
{
    return x < 4 ? 0 : floorlog4(x - 3) + 1;
}

/**
 * @brief Compile time floorlog8
 * 
 * @param x Input
 * @return constexpr unsigned floor(log8(x))
 */
constexpr unsigned floorlog8(unsigned x)
{
    return x < 8 ? 0 : 1+floorlog8(x >> 3);
}

/**
 * @brief Compile time ceil log8
 * 
 * @param x Input
 * @return constexpr unsigned ceil(log8(x))
 */
constexpr unsigned ceillog8(unsigned x)
{
    return x < 8 ? 0 : floorlog8(x - 7) + 1;
}

/**
 * @brief Compile time ceil log radix R
 * 
 * @param x input value
 * @param R radix
 * @return constexpr unsigned ceil(logR(x))
 */
constexpr unsigned ceillogR(unsigned x, unsigned R){
    if(R == 2){
        return ceillog2(x);
    }
    if(R == 4){
        return ceillog4(x);
    }
    if(R == 8){
        return ceillog8(x);
    }
    return 0;
}

#endif // CONST_LOGS_HPP_