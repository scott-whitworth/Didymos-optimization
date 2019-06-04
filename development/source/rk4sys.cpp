
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

            if(n<4)
            {
                std::cout<<"Y state "<<n<<std::endl<<std::endl;
                std::cout<<"r= "<<y[n].r<<std::endl;
                std::cout<<"theta= "<<y[n].theta<<std::endl;
                std::cout<<"vr= "<<y[n].vr<<std::endl;
                std::cout<<"vtheta= "<<y[n].vtheta<<std::endl<<std::endl;
            }
            k1 = calc_k(y[n],stepSize/2);
            ymid = calc_ymid(y[n],stepSize/2,k1);
             if(n<4)
            {
                std::cout<<"K1 - loop "<<n<<std::endl<<std::endl;
                std::cout<<"rmid= "<<ymid.r<<std::endl;
                std::cout<<"thetamid= "<<ymid.theta<<std::endl;
                std::cout<<"vrmid= "<<ymid.vr<<std::endl;
                std::cout<<"vthetamid= "<<ymid.vtheta<<std::endl;
                std::cout<<"K1-r = "<<k1.r<<std::endl;
                std::cout<<"K1-theta = "<<k1.theta<<std::endl<<std::endl;
            }
            k2 = calc_k(ymid,stepSize/2);
            ymid = calc_ymid(ymid,stepSize/2,k2);
            if(n<4)
            {
                std::cout<<"K2 - loop "<<n<<std::endl<<std::endl;
                std::cout<<"rmid= "<<ymid.r<<std::endl;
                std::cout<<"thetamid= "<<ymid.theta<<std::endl;
                std::cout<<"vrmid= "<<ymid.vr<<std::endl;
                std::cout<<"vthetamid= "<<ymid.vtheta<<std::endl;
                std::cout<<"K2-r = "<<k2.r<<std::endl;
                std::cout<<"K2-theta = "<<k2.theta<<std::endl<<std::endl;
            }
            k3 = calc_k(ymid,stepSize);
            ymid = calc_ymid(ymid,stepSize,k3);
             if(n<4)
            {
                std::cout<<"K3 - loop "<<n<<std::endl<<std::endl;
                std::cout<<"rmid= "<<ymid.r<<std::endl;
                std::cout<<"thetamid= "<<ymid.theta<<std::endl;
                std::cout<<"vrmid= "<<ymid.vr<<std::endl;
                std::cout<<"vthetamid= "<<ymid.vtheta<<std::endl;
                std::cout<<"K3-r = "<<k3.r<<std::endl;
                std::cout<<"K3-theta = "<<k3.theta<<std::endl<<std::endl;
            }
		    k4 = calc_k(ymid, stepSize);

            // Add weighted slopes
			elements phi = (k1 + (k2 + k3) * 2 + k4) / 6; // calculate phi for each element
            y[n+1] = y[n] + phi * stepSize;
        

    }

    return y;
}