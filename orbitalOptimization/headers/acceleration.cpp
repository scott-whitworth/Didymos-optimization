#include "acceleration.h"
#include "constants.h"
#include "thruster.h"
#include <iostream>

template <class T> T calc_accel(const T & r, thruster<T> & thrusterType,T & massFuel, const T & deltaT, const bool & thrusting){
    
    if(!thrusting){
        return 0;
    }

    if(wetMass - massFuel < dryMass){
        return 0;
    }
    
    T Pin;
    T Pthrust;
    T thrust;
    T m_spaceCraft;
    T accel;

    Pin = thrusterType.P0/pow(r,2);

    Pthrust = thrusterType.calc_eff(Pin)*Pin;

    thrusterType.calc_m_Dot(Pin);

    thrust = sqrt(2*Pthrust*thrusterType.m_Dot);

    massFuel = massFuel + thrusterType.m_Dot*deltaT;

    m_spaceCraft = wetMass - massFuel;
    
    accel = thrust/m_spaceCraft;

    accel = accel/AU;

return accel;
}

