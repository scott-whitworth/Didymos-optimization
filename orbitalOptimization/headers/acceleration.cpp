#include "acceleration.h" 
#include "constants.h" // used for wetMass
#include "thruster.h" // used for P0, calc_eff(), and calc_m_Dot()
#include <iostream> // used for cout

template <class T> T calc_accel(const T & radius, thruster<T> & thrusterType, T & massFuel, const T & deltaT, const T & thrusting, const T & wetMass){
    
    // Ensures there is a limit to how much fuel is used in a given trip.
    // wetMass is a optimized variable and dryMass is defined in constants.h 
    if(wetMass - massFuel < dryMass){
        return 0;
    }

    // Thrusting is evaluated in calcFourier.cpp within calc_coast().
    // When thrusting is equal to zero, calc_accel() will not be evaluated.
    if(thrusting == 0){
        return 0;
    }

    //Defining variables for calc_accel().
    T Pin;
    T Pthrust;
    T thrust;
    T m_spaceCraft;
    T accel;

    // Power going into the spacecraft as a function of the radius of the spacecraft from the sun.
    Pin = thrusterType.P0/pow(radius,2);

    //If the spacecraft is closer to the sun than the earth, the power in can not be greater than the power measured on earth.
    if(radius<1){
        Pin = thrusterType.P0/1;
    }

    //The thrust power of the spacecraft is dependent upon the efficiency (calculated in thruster.cpp) and the power (in).
    Pthrust = thrusterType.calc_eff(Pin)*Pin;

    //Calling calc_m_Dot (defined in thruster.cpp) which returns the fuel flow rate (mDot) for a given power (in).
    thrusterType.calc_m_Dot(Pin);

    //Thrust is calculated by power (thrust) and mDot.
    thrust = sqrt(2*Pthrust*thrusterType.m_Dot);

    //Calculates the amount of fuel used throughout the duration of the trip.
    massFuel = massFuel + thrusterType.m_Dot*deltaT;

    //Calculates the current mass of the spacecraft given the amount of fuel used subtracted from the wetMass(defined in constants.h).
    m_spaceCraft = wetMass - massFuel;
    
    //Acceleration of the spacecraft due to thrusting calculated by thrust divided by the mass of the spacecraft.
    // AU converts the acceleration from m/s^2 to au/s^2.
    accel = thrust/(AU*m_spaceCraft);

return accel;
}

