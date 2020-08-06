#ifndef thruster_h
#define thruster_h

// File path for all data used in this file: \\fs1\phys\sankaranResearch\2019-Lauren-Mateo\Thruster

// Sets starting values as given in the 2017 data excel sheet
template <class T> struct thruster {
    T m_Dot; // Fuel flow rate
    int type; // Type of thruster
    T P0; // Initial power in
    T m_Dot0; // Initial m_Dot

    thruster<T>(int type);

    // Overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const thruster<T> & e); 

// Evaluates the spacecraft's effciency for a certian iteration based off of a best fit line from file above of eff vs. power in (W).
// Parameters:
//         Pin: Given a certian value of power in to the spacecraft, a specific efficiency will be calculated for an iteration.
// Output: spacecraft's effciency for a certian iteration
    T calc_eff(const T & Pin);

// Evaluates the spacecraft's fuel flow rate (mDot) for a certian iteration based off of "if statements".
// Parameters:
//         Pin: Given a certian value of power in to the spacecraft, a certian fuel flow rate will be set for an iteration.
// Output: spacecraft's mDot for a certian iteration
    void calc_m_Dot(const T & Pin);

    private:
    T NEXTP0 =7330; // Initial power (W)
    T NEXTm_Dot0 = 5.73E-06; // Initial fuel flow rate (kg/s)
};

#include "thruster.cpp"
#endif