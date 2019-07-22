#include "optimizationCoeff.h"
#include <iomanip> // setprecision(int)

OptimizationCoeff::OptimizationCoeff(){
    for(int i = 0; i < gammaSize; i++){
        gamma[i] = 1;
    }
    for(int i = 0; i < tauSize; i++){
        tau[i] = 2;
    }
    for(int i = 0; i < coastSize; i++){
        coast[i] = 3;
    }

    coastThreshold = 15;
}

std::ostream & operator<<(std::ostream & Str, const OptimizationCoeff & e){
    Str << std::fixed;
    Str << std::setprecision(4); // number of decimals output into text file
    Str << "Gamma [";
    for(int i=0;i<e.gammaSize;i++){
        Str << e.gamma[i];
        if(i < e.gammaSize-1){
            Str << ", ";
        }
    }
    Str << "]" << std::endl;

    Str << "tau [";
    for(int i=0;i<e.tauSize;i++){
        Str << e.tau[i];
        if(i < e.tauSize-1){
            Str << ", ";
        }
    }
    Str << "]" << std::endl;

    Str << "Coast [";
    for(int i=0;i<e.coastSize;i++){
        Str << e.coast[i];
        if(i < e.coastSize-1){
            Str << ", ";
        }
    }
    Str << "]" << std::endl;
    Str << "Coasting Threshold: " << e.coastThreshold << std::endl;

    return Str;
}