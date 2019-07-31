//Didymos-Optimization_Project:
//Last Editor: Ben and Mateo
//Tasks Completed: 
    //Changed the magic numbers to defined constants
#ifndef coefficients_cpp
#define coefficients_cpp

#include "coefficients.h"
#include <iomanip> // setprecision(int)

template <class T> std::ostream & operator<<(std::ostream & Str, const coefficients<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << "Gamma: \t";
    for(int i = 0; i < e.gammaSize; i++){
        Str << e.gamma[i];
        if(!(i == e.gammaSize-1) ){
            Str << ",\t";
        }
    }
    Str << endl << "Tau: \t";
    for(int i = 0; i < e.tauSize; i++){
        Str << e.tau[i];
        if(!(i == e.tauSize-1) ){
            Str << ",\t";
        }
    }
    Str << endl << "Coasting: \t";
    for(int i = 0; i < e.coastSize; i++){
        Str << e.coast[i];
        if(!(i == e.coastSize-1) ){
            Str << ",\t";
        }
    }
    Str << endl;
    // TODO: fix coasting threshold to be a constant
    Str << "Coasting Threshold: " << e.coastThreshold << endl;
    return Str;
}

#endif