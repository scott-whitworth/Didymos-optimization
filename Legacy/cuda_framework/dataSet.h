#ifndef DATASET_H
#define DATASET_H

#include <iostream>
#include "d_trip.h"
#include "optimizationCoeff.h"

class DataSet{

    public:
        D_Trip pos_vel;
        OptimizationCoeff coeff;

        //Constructors
        DataSet();
        DataSet(D_Trip in_pos_vel, OptimizationCoeff in_coeff);

        friend std::ostream & operator<<(std::ostream & Str, const DataSet & e); 
};

#endif