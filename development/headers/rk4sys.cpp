//#include "rk4sys.h"
#include <iostream>
#include "elements.h"

template <class T> elements<T>* rk4sys(T timeInitial, T timeFinal, elements<T> y0, T stepSize, int numSteps, elements<T> *y){
    // Set the first element of the solution vector to the initial conditions
    y[0] = y0;
    // k variables for Runge-Kutta calculation of y[n+1]
    elements<T> k1, k2, k3, k4;

    for(int n=0;n<numSteps-1;n++) // iterate over all time steps 
    {
        // If we required the time
        //      t = the number of iterations multiplied by the time interval between points
       
       // Runge-Kutta algorithm
        k1 = calc_k(stepSize, y[n]); //h*y[n]
    
        k2 = calc_k(stepSize, y[n]+k1/2); //h*(y[n]+k1/2)

        k3 = calc_k(stepSize, y[n]+k2/2); //h*(y[n]+k2/2)

        k4 = calc_k(stepSize, y[n]+k3); //h*(y[n]+k3)

        // Add weighted slopes (k elements)
        y[n+1] = y[n] + (k1 + (k2 + k3) * 2 + k4) / 6; // calculates the y[n] for the next round of calculations
    }
    return y;
}