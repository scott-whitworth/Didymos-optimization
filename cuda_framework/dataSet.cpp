#include "dataSet.h"

//Constructors
DataSet::DataSet() : pos_vel(), coeff(){

}

DataSet::DataSet(D_Trip in_pos_vel, OptimizationCoeff in_coeff) : pos_vel(in_pos_vel), coeff(in_coeff){

}

std::ostream & operator<<(std::ostream & Str, const DataSet & e){
    Str << e.pos_vel << std::endl;
    Str << e.coeff;
    return Str;
}