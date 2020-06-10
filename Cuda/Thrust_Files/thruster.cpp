//Didymos-Optimization_Project:
//Last Editor: Ben and Lauren
//Tasks Completed: 
      //Created thruster file which sets the thruster being used in the optimization.
      //NEXT:
        //Inserted a function which relates effciency to the power going in to the spacecraft based off of an excel best fit line.
        //Created a variable mDot by adding if statements which relate the power going in to the spacecraft to the amount of fuel consumption for a given iteration.

#include "thruster.h"
#include <iomanip> // used for setPrecision()

template <class T>
thruster<T>::thruster(int newType) {
    // setting values (defined in header) for when type 1 is called (NEXT)
    if (newType == 1) {
        m_Dot = NEXTm_Dot0;
        P0 = NEXTP0;
    }
    type = newType;
}

template <class T> T thruster<T>::calc_eff(const T & Pin) {
   // Data interpolation for thruster type NEXT
   if (type == 1) {
    return  -1.328086e-23*pow(Pin,6) + 6.207694e-19*pow(Pin,5) - 9.991813e-15*pow(Pin,4) +  7.701266e-11*pow(Pin,3) - 3.136031e-07*pow(Pin,2) +  6.805225e-04*Pin;       // Polynomial fit
   }
   return 0;
}

template <class T> void thruster<T>::calc_m_Dot(const T & Pin) {
    if (type == 1) {
        if (Pin < 2550) {
            this->m_Dot = 1.99E-06;
        }
        else if (Pin < 4500) {
            this->m_Dot = 4.44E-06;
        }
        else {
            this->m_Dot = 5.73E-06;
        }
    }
}


template <class T> std::ostream & operator<<(std::ostream & Str, const thruster<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.gamma << "\t" << e.tau << "\n";
    return Str;
}