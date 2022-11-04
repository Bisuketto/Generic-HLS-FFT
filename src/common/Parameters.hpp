/**
 * @file Parameters.hpp
 * @author Bordeaux INP, Bertrand LE GAL (bertrand.legal@ims-bordeaux.fr) [http://legal.vvv.enseirb-matmeca.fr]
 * @brief This file is part of a LDPC library for realtime LDPC decoding on processor cores.
 * @version 0.0.0
 * @date 2020
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#ifndef CLASS_Parameters
#define CLASS_Parameters

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <map>

using namespace std;

class Parameters{
	bool locked;
	map<string, string> liste;

public:
    Parameters();
	~Parameters();

protected:
	string param(const string p);

public:
	void set(const string p, const string   v);
	void set(const string p, const int      v);
	void set(const string p, const uint32_t v);
	void set(const string p, const double   v);
	void set(const string p, const float    v);
//	void set(const string p, const bool     v);

    bool exist(const string p);

    string toString(const string p);
    int    toInt   (const string p);
    long   toLong  (const string p);
    float  toFloat (const string p);
    double toDouble(const string p);
    bool   toBool  (const string p);
    string toBoolString(const string p);

    };

#endif // CLASS_CTools
