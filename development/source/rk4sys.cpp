
#include <vector>
#include "rk4sys.h"
#include <iostream>

/* fourth-order RUnge-Kutta for a system of ODEs
-integrates a system of ODEs with fourth-order RK method

Input:
double timeInitial and double timeFinal; initial and final times for the computation,
elements y0: initial values of dependent variables
stepSize: the difference between times
output:
returns: y - solutions of dependent variables
*/ 
elements* rk4sys(double timeInitial, double timeFinal, elements y0, double stepSize){
    // Define the max number of iterations
    int nMax = (int) (((timeFinal-timeInitial)/stepSize)+0.5);

    // How to allocate memory in C
    // elements* y;
    //  y = (elements *)malloc(sizeOf(elements)*nMax);
    elements* y;
    y = new elements[nMax];
    // Set the first element of the solution vector to the initial conditions
    y[0] = y0;
    
    for(int n=0;n<nMax-1;n++)
    {
        // If we required the time
        // time = stepSize*n

        // Variables for Runge-Kutta
        elements k1, k2, k3, k4, ymid;
       
       // Runge-Kutta algorithm

            k1 = calc_k(stepSize, y[n]);
        
            k2 = calc_k(stepSize, y[n]+k1/2);

            k3 = calc_k(stepSize, y[n]+k2/2);

            k4 = calc_k(stepSize, y[n]+k3);

            // Add weighted slopes
			elements phi = (k1 + (k2 + k3) * 2 + k4) / 6; // calculate phi for each element
            y[n+1] = y[n] + phi;
    }
    return y;
}