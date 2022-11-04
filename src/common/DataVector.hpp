/**
 * @file DataVector.hpp
 * @author Bordeaux INP, Bertrand LE GAL (bertrand.legal@ims-bordeaux.fr) [http://legal.vvv.enseirb-matmeca.fr]
 * @brief This file is part of a LDPC library for realtime LDPC decoding on processor cores.
 * @version 0.0.0
 * @date 2020
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef _DataVector_
#define _DataVector_

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <fstream>
#include <map>
#include <sstream>
#include <cmath>
#include <complex>

using namespace std;

class DataVector{
public:
    std::vector<float> I;
    std::vector<float> Q;

    DataVector(const DataVector& in);

    DataVector(const std::string filename);

    bool save(const std::string filename);

//    bool copyInputsToOutputs();
};

#endif
