#ifndef thruster_h
#define thruster_h
#include "constants.h"



template <class T> struct thruster {
    T m_Dot;
    int type;
    T P0;
    T m_Dot0;

    //thruster<T>(T m_Dot0, T eff0, T P00);
    thruster<T>(int type);

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const thruster<T> & e); 
    T calc_eff(const T & Pin);
    void calc_m_Dot(const T & Pin);
    private:
    T NEXTP0 =7330; // initial power (W)
    T NEXTm_Dot0 = 5.73E-06;
};

#include "thruster.cpp"
#endif