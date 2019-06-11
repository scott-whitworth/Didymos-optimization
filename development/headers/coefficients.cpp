#include "coefficients.h"
#include <iomanip> // setprecision(int)


template <class T> std::ostream & operator<<(std::ostream & Str, coefficients<T> const & e) {
    Str << std::fixed;
    Str << std::setprecision(16); // number of decimals output into text file
    Str << e.gamma << "\t" << e.tau << "\n";
    return Str;
}