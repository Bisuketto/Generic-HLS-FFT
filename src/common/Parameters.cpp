/**
 * @file Parameters.cpp
 * @author Bordeaux INP, Bertrand LE GAL (bertrand.legal@ims-bordeaux.fr) [http://legal.vvv.enseirb-matmeca.fr]
 * @brief This file is part of a LDPC library for realtime LDPC decoding on processor cores.
 * @version 0.0.0
 * @date 2020
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include "Parameters.hpp"

Parameters::Parameters() : liste(){
	locked = false;
	liste["snr_min"         ] = "0.50";
	liste["snr_max"         ] = "4.01";
	liste["snr_pas"         ] = "0.25";
	liste["fe_limit"        ] = "100";
	liste["channel_type"    ] = "3";
//	liste["poinconnage"     ] = "0";
	liste["norm_channel"    ] = "false";
	liste["real_encoder"    ] = "false";
	liste["qpsk_channel"    ] = "false";
	liste["Es_N0"           ] = "false";
	liste["worst_case_fer"  ] = "false";
	liste["show_llr_histo"  ] = "false";
	liste["refresh"         ] = "1";
}

Parameters::~Parameters(){

}

string Parameters::param(const string p){
	string r = liste[p];
    if( r == "" )
    {
        std::cout << "(EE) Une erreur a été detectée, le parametre (" << p << ") n'existe pas !" << std::endl;
        exit( EXIT_FAILURE );
    }
	assert( r != "" );
	return r;
}

void Parameters::set(const string p, const string v){
	assert( locked == false );
	liste[p] = v;
}

void Parameters::set(const string p, const int v){
	string w = to_string(v);
	set(p, w);
//	cout << "SET Parameters[" << p << "] = " << w << " (" << liste[p] << ")" << endl;
}

void Parameters::set(const string p, const uint32_t v){
	set(p, to_string(v));
}

void Parameters::set(const string p, const double v){
	set(p, to_string(v));
}

void Parameters::set(const string p, const float v){
	set(p, to_string(v));
}
/*
void Parameters::set(const string p, const bool v)
{
	string w = v ? "true" : "false";
	set(p, w);
	
	if (v == true) {
		set(p, 1);
	} else {
		set(p, 0);
	}

*/
bool Parameters::exist(const string p){
	string r = liste[p];
	return (r != "");
}

string Parameters::toString(const string p){
	return param(p);
}

int Parameters::toInt(const string p){
	string w = param(p);
	return stoi( w );
}

long Parameters::toLong(const string p){
	return stol( param(p) );
}

float Parameters::toFloat(const string p){
	return stof( param(p) );
}

double Parameters::toDouble(const string p){
	return stod( param(p) );
}

bool Parameters::toBool(const string p){
	return toInt( p );
	/*
	string v = param(p);
	if( v ==  "true" ) return true;
	if( v ==     "1" ) return true;
	if( v == "false" ) return false;
	if( v ==     "0" ) return false;
	assert( false );
	return false;
*/
}

std::string Parameters::toBoolString(const string p){
    if( toInt( p ) == 1 )
        return "on";
    else
        return "off";
}
