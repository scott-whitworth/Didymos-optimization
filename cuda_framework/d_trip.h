#ifndef D_TRIP_H
#define D_TRIP_H

#include <iostream>


class D_Trip{

    public:
        //Position Data
        double r;
        double theta;
        double z;

        //Velocity Data
        double v_r;
        double v_theta;
        double v_z;

        //Constructors
        D_Trip();
        D_Trip(double in_r, double in_theta, double in_z, double in_v_r, double in_v_theta, double in_v_z);
        D_Trip(double in_r, double in_theta, double in_z);

        friend std::ostream & operator<<(std::ostream & Str, const D_Trip & e); 
};


#endif