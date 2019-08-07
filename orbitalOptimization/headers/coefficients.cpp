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
    Str << e.gamma << "\t" << e.tau << "\n";
    return Str;
}

template <class T> void initCoefficient(double x[], coefficients<T> & coeff)
{
    // in-plane angle
    for (int i=0;i<coeff.gammaSize;i++)
    {
        coeff.gamma[i]=x[i+GAMMA_OFFSET];
    }
    // out-of-plane angle
    for (int i=0;i<coeff.tauSize;i++)
    {
        coeff.tau[i]=x[i+TAU_OFFSET];
    }

    // setup of coast determination calculations based off of optimized coefficients
    for (int i=0;i<coeff.coastSize;i++)
    {
        coeff.coast[i]=x[i+COAST_OFFSET];
    }
}

#endif