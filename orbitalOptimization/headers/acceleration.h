#ifndef acceleration_h
#define acceleration_h

#include "thruster.h"

template <class T> T calc_accel(const T & r, const thruster<T> & thrusterType, const T & massFuel, const T & deltaT, const T & thrusting, const T & dryMass);

#include "acceleration.cpp"
#endif