#include "acceleration.h"
#include "constants.h"
#include <iostream>

template <class T> T calc_accel(const T & r, const thruster<T> & thrusterType, T & massFuel, const T & deltaT){

    T Pin;
    T Pthrust;
    T thrust;
    T m_spaceCraft;
    T accel;

    Pin = thrusterType.P0/pow(r,2);

    Pthrust = thrusterType.eff*Pin;

    thrust = sqrt(2*Pthrust*thrusterType.m_Dot);

    massFuel = massFuel + thrusterType.m_Dot*deltaT;

    m_spaceCraft = wetMass - massFuel;
    
    accel = thrust/m_spaceCraft;

    accel = accel/AU;

return accel;
}
