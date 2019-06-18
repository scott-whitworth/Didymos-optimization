#ifndef thruster_h
#define thruster_h



template <class T> struct thruster {

    T m_Dot;
    T eff;
    T P0;

    thruster<T>(T m_Dot0, T eff0, T P00);

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const thruster<T> & e); 
};

#include "thruster.cpp"
#endif