
#include "acceleration.h" 
#include <iostream> // used for cout

template <class T> __host__ __device__ T calc_accel(const T & radius, const T & z, thruster<T> & thrusterType, T & massExpelled, const T & deltaT, const bool & thrusting, const T & wetMass, const cudaConstants* cConstants) {

    // If all of the fuel has been expelled, then no more thrust can be applied
    if (wetMass - massExpelled <= cConstants->dry_mass) {
        return 0;
    }

    // Thrusting is evaluated in calcFourier.cpp within calc_coast().
    // When thrusting is equal to zero, calc_accel() will not be evaluated.
    if (!thrusting) {
        return 0;
    }

    // Defining variables for calc_accel().
    T Pin; // Power input
    T Pthrust; // Thrust power
    T thrust;

    // Power going into the spacecraft as a function of the radius of the spacecraft from the sun (r is non-dimensionalized by dividing by 1 AU).
    T sepparation = sqrt( pow(radius, 2) + pow(z, 2) );
    Pin = thrusterType.P0/sepparation; 

    //If the spacecraft is closer to the sun than the earth, the power in can not be greater than the power experimentally measured on earth.
    //This creates a "sphere" around the sun to ensure the power does not exceed the tested limit.
    if (sepparation <= 1) {
        Pin = thrusterType.P0; // It is devided by 1 astronomical unit to normalize it P0/(1 AU)
    }

    // The thrust power of the spacecraft is dependent upon the efficiency (calculated in thruster.cpp) and the power (in).
    Pthrust = thrusterType.calc_eff(Pin)*Pin; 

    // Update thrusterType's current m_Dot based on power input
    thrusterType.calc_m_Dot(Pin);

    // Thrust is calculated by power (thrust) and mDot.
    thrust = sqrt(2 * Pthrust * thrusterType.m_Dot); 

    // Calculates the amount of fuel used throughout the duration of the trip.
    massExpelled += thrusterType.m_Dot * deltaT;
    
    // the current mass of the spacecraft is equal to the fuel used minus the wetMass of the spacecraft
    // Acceleration of the spacecraft due to thrusting calculated by thrust divided by the mass of the spacecraft.
    // AU converts the acceleration from m/s^2 to au/s^2.
    return thrust/(AU*(wetMass - massExpelled));
}
