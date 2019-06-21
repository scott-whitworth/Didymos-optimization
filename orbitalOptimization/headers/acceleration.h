#ifndef acceleration_h
#define acceleration_h

#include "thruster.h"

template <class T> T calc_accel(const T & r, thruster<T> & thrusterType, T & massFuel, const T & deltaT, const bool & thrusting);

#include "acceleration.cpp"
#endif