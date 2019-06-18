#include "thruster.h"
#include <iomanip> // setprecision(int)


// sets starting values as given
template <class T>
thruster<T>::thruster(T m_Dot0, T eff0, T P00){
    m_Dot = m_Dot0;
    eff = eff0;
    P0 = P00;
}

template <class T> std::ostream & operator<<(std::ostream & Str, const thruster<T> & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.gamma << "\t" << e.tau << "\n";
    return Str;
}