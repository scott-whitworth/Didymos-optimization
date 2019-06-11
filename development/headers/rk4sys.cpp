#include "rk4sys.h"
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> elements<T>* rk4sys( const T & timeInitial, const T & timeFinal, T *times, const  elements<T> & y0, T & stepSize, elements<T> *y, const T & absTol, coefficients<T> & coeff, const T & accel, T *gamma, T *tau){
    // Set the first element of the solution vector to the initial conditions
    y[0] = y0;
    times[0]=timeInitial;
    // k variables for Runge-Kutta calculation of y[n+1]
    elements<T> k1, k2, k3, k4, k5, k6, k7;

    T t = timeInitial; // setting time equal to the start time
    int n=0; // setting the initial iteration number equal to 0

    //TODO: SC: This while loop (and the function in general) should be functionalized. This whole function has gotten to the point where it is sufficiently complex to cause issues moving forward

    while(t<=timeFinal) // iterate until time is equal to the stop time
    {
//      array of gamma for binary output
        gamma[n] =calc_gamma(coeff,t, timeFinal);
//      array of tau for binary output
        tau[n] =calc_tau(coeff,t, timeFinal);  
        // array of time output as t         
        t += stepSize;
        // Time of iteration is set to the previous time plus the step size used within that iteration
        times[n+1]=t;

        // Runge-Kutta algorithm       
        // k1 = h*f(t, y[n])
        k1 = calc_k(stepSize, y[n], coeff, accel, t, timeFinal);        
        // k2 = h*f(t+1/5, y[n]+k1*1/5)
        k2 = calc_k(stepSize, y[n]+k1*1/5,coeff, accel, t+1/5*stepSize, timeFinal);   
        // k3 = h*f(t+3/10, y[n]+k1*3/40+k2*9/40)
        k3 = calc_k(stepSize, y[n]+k1*3/40+k2*9/40,coeff, accel, t+3/10*stepSize, timeFinal);   
        // k4 = h*f(t+4/5, y[n]+k1*44/45+k2*-56/15+k3*32/9)
        k4 = calc_k(stepSize,y[n]+k1*44/45+k2*-56/15+k3*32/9,coeff, accel, t+4/5*stepSize, timeFinal);    
        // k5 = h*f(t+8/9, y[n]+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729)
        k5 = calc_k(stepSize, y[n]+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729,coeff, accel, t+8/9*stepSize, timeFinal);        
        // k6 = h*f(t, y[n]+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656)
        k6 = calc_k(stepSize, y[n]+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656,coeff, accel, t+stepSize, timeFinal);        
        // k7 = h*f(t, y[n]+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84)
        k7 = calc_k(stepSize,y[n]+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84,coeff, accel, t+stepSize, timeFinal);  

        // Previous value 
        // v = y[n] + 5179/57600*k1 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7
        elements<T> v = y[n] + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;     

        // Current value
        // u = y[n] + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
        elements<T> u = y[n] + k1*(35./384) + k3*(500./1113) + k4*125./192 - k5*2187/6784 + k6*11/84;  

        // Alter the step size for the next iteration
        T s = calc_scalingFactor(v,u-v,absTol,stepSize);
        stepSize = s*stepSize;

        //TODO: SC: You take a slightly risky move here in changing stepSize. I know it makes sense in the context of the function, but it is also a variable you pass in.
        //          This discrepancy might cause issues moving forward

        // The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/2){
            stepSize=(timeFinal-timeInitial)/2;
        } else if (stepSize<(timeFinal-timeInitial)/1000){
            stepSize=(timeFinal-timeInitial)/1000;
        }

        //TODO: SC: This comment is not what it is doing:
        // Calculates the y[n] for the next round of calculations
        y[n+1] = u;        
        n++;
    }//end of while 

    return y; //TODO: SC: Why are we returning y?
}

template <class T> T calc_scalingFactor(const elements<T> & previous, const elements<T> & difference, const T & absTol, const T & stepSize)
{
    //TODO: SC: normTotError might mean Normilize Total Error, but it should be documented. Also what does scale refer to? Why do we need to variables?
    T normTotError, scale;

        // relative error (unitless) 
        elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
        difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);

    //TODO: SC: This comment is not complete, and not accurate. You are taking the square root of the sum of squares of the errors
    //TODO: SC: using (T) as a casting method leads to some issues. The better way is static_cast<T>(var). I am also not totally sure you need to cast anything here. The pmError elements will need a T pow function call
    // Sum the error of the 6 element to determine the scale for the time step of the next iteration
    normTotError = pow(pow(pmError.r,(T)2) + pow(pmError.theta,(T)2) + pow(pmError.z,(T)2) + pow(pmError.vr,(T)2) + pow(pmError.vtheta,(T)2) + pow(pmError.vz,(T)2),(T)1/2);
    scale = pow((absTol/normTotError),(T)1/5);

        return scale;   
}