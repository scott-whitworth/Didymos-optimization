//Didymos-Optimization_Project:

#include "thruster.h"
#include <iomanip> // used for setPrecision()

template <class T>
thruster<T>::thruster(const cudaConstants* gConfig) {

    // If no thruster, set values to 0
    if (gConfig->thruster_type == THRUST_TYPE::NO_THRUST) {
        m_Dot = P0 = 0;
    }

    // setting values (defined in header) for when type 1 is called (NEXT)
    else if (gConfig->thruster_type == THRUST_TYPE::NEXT_C) {
        m_Dot = m_Dot0 = NEXTm_Dot0;
        P0 = NEXTP0;
    }
    
    type = gConfig->thruster_type;
    coastThreshold = gConfig->coast_threshold;
}

template <class T> T thruster<T>::calc_eff(const T & Pin) {
    // Data interpolation for thruster type NEXT
    if (type == THRUST_TYPE::NEXT_C) {
        return  -1.328086e-23*pow(Pin,6) + 6.207694e-19*pow(Pin,5) - 9.991813e-15*pow(Pin,4) +  7.701266e-11*pow(Pin,3) - 3.136031e-07*pow(Pin,2) +  6.805225e-04*Pin;       // Polynomial fit
    }
    return 0;
}

template <class T> void thruster<T>::calc_m_Dot(const T & Pin) {
    if (type == THRUST_TYPE::NEXT_C) {
        if (Pin < 2550) {
            m_Dot = 1.99E-06;
        }
        else if (Pin < 4500) {
            m_Dot = 4.44E-06;
        }
        else {
            m_Dot = NEXTm_Dot0;
        }
    }
}


template <class T> std::ostream & operator<<(std::ostream & Str, const thruster<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.m_Dot << "\t" << e.P0 << "\n";
    return Str;
}
