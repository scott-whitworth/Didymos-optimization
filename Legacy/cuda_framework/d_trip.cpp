#include "d_trip.h"
#include <iomanip> // setprecision(int)

D_Trip::D_Trip() : D_Trip(0.0, 0.0, 0.0, 0.0, 0.0, 0.0){

}

D_Trip::D_Trip(double in_r, double in_theta, double in_z, double in_v_r, double in_v_theta, double in_v_z){
    r = in_r;
    theta = in_theta;
    z = in_z;

    v_r = in_v_r;
    v_theta = in_v_theta;
    v_z = in_v_z;
}

D_Trip::D_Trip(double in_r, double in_theta, double in_z) : D_Trip(in_r, in_theta, in_z, 0.0, 0.0, 0.0) {
    
}

std::ostream & operator<<(std::ostream & Str, const D_Trip & e){
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.r << "\t" << e.theta << "\t" << e.z << "\t" << e.v_r << "\t" << e.v_theta << "\t" << e.v_z;
    return Str;
}