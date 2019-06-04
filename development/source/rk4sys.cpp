
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
    int nMax = (timeFinal-timeInitial)/stepSize;

    // How to allocate memory in C
    // elements* y;
    //  y = (elements *)malloc(sizeOf(elements)*nMax);
    elements* y;
    y = new elements[nMax];
    // Set the first element of the solution vector to the initial conditions
    y[0] = y0;
    
    for(int n=0;n<=nMax+1;n++)
    {
        // If we required the time
        // time = stepSize*n
        
        // Variables for Runge-Kutta
        elements k1, k2, k3, k4, ymid;
       
       // Runge-Kutta algorithm

            k1 = calc_k(y[n],stepSize/2);
            ymid = calc_ymid(ymid,stepSize,k1);
            k2 = calc_k(y0,stepSize/2);
            ymid = calc_ymid(ymid,stepSize,k2);
            k3 = calc_k(ymid,stepSize);
            ymid = calc_ymid(ymid,stepSize,k3);
		    k4 = calc_k(ymid, stepSize);

            // Add weighted slopes
			elements phi = (k1 + (k2 + k3) * 2 + k4) / 6; // calculate phi for each element
            y[n+1] = y[n] + phi * stepSize;
        std::cout<<y[n].r<<std::endl;
        std::cout<<k1.r<<std::endl;

    }

    return y;
}