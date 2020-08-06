#ifndef OPTIMIZATIONCOEFF_H
#define OPTIMIZATIONCOEFF_H

#include <iostream>

class OptimizationCoeff{

    public:
        const static int gammaSize = 9;
        const static int tauSize = 5;
        const static int coastSize = 3;

        double gamma[gammaSize];
        double tau[tauSize];
        double coast[coastSize];

        double coastThreshold;

        OptimizationCoeff();

        friend std::ostream & operator<<(std::ostream & Str, const OptimizationCoeff & e); 

};

#endif