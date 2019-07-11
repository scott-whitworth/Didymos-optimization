#include "orbitalMotion.h"
#include "elements.h"

#include <iostream>

int main()
{
    elements<double> result;
    result = earthInitial(2.99017e+007);
    std::cout<<result<<std::endl;
}