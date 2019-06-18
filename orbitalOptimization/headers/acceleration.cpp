#include "acceleration.h"
#include "constants.h"

template <class T> T calc_accel(const T & r, const thruster<T> & thrusterType, const T & massFuel, const T & deltaT){

    massFuel = massFuel + thrusterType.m_dot*deltaT;

    Pin = thrusterType.P0/r^2;

    Pthrust = thrusterType.eff*Pin;

    thrust = sqrt(2*Pthrust*thrusterType.mdot)

    m_spaceCraft = wetMass - massFuel;
    
    accel = thrust/m_spaceCraft;

}
