#include "rk4sys.h"
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y, const T & absTol, coefficients<T> coeff, const T & accel, T *gamma,  T *tau){
<<<<<<< HEAD
    
=======
        // Set the first element of the solution vector to the initial conditions
        y[0] = y0;
        times[0]=timeInitial;
>>>>>>> b138334863b971538f35f81cb1ae10de86533b03

        T curTime = timeInitial; // setting time equal to the start time
        int iteration=0; // setting the initial iteration number equal to 0

    //TODO: SC: This while loop (and the function in general) should be functionalized. This whole function has gotten to the point where it is sufficiently complex to cause issues moving forward

        while(curTime<=timeFinal) // iterate until the current time reaches the end time
        {
                rkStep(&curTime, timeInitial, timeFinal, times, stepSize, y, absTol, coeff, accel, gamma, tau, &iteration); //take a step in the runge-kutta algorithm
        }
}

template <class T> void rkStep(T *curTime, const T & timeInitial, const T & timeFinal, T *times, T stepSize, elements<T> *y, const T & absTol, coefficients<T> coeff, const T & accel, T *gamma, T *tau, int *iteration){
        elements<T> *v; // used to store the previous and current values in rkCalc()
        elements<T> *u;

        // array of gamma for binary file output
        gamma[*iteration] = calc_gamma(coeff, *curTime, timeFinal);
        // array of tau for binary file output
        tau[*iteration] = calc_tau(coeff, *curTime, timeFinal);  

        // array of time output as t         
        *curTime += stepSize;
        // Time of iteration is set to the previous time plus the step size used within that iteration
        times[*iteration+1]=*curTime;

        rkCalc(curTime, timeFinal, stepSize, y, coeff, accel, v, u, *iteration); //calculate new elements values for the current time

        // Alter the step size for the next iteration
        T s = calc_scalingFactor(*v,*u-*v,absTol,stepSize);
        stepSize = s*stepSize;

        //TODO: SC: You take a slightly risky move here in changing stepSize. I know it makes sense in the context of the function, but it is also a variable you pass in.
        // This discrepancy might cause issues moving forward

        // The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/2)
                stepSize=(timeFinal-timeInitial)/2;
        else if (stepSize<(timeFinal-timeInitial)/1000)
                stepSize=(timeFinal-timeInitial)/1000;

        // store the new set of elements
        y[*iteration+1] = *u;        
        *iteration++;
}

template <class T> void rkCalc(T *curTime, const T & timeFinal, T stepSize, elements<T> *y, coefficients<T> coeff, const T & accel, elements<T> *v, elements<T> *u, int iteration){
        // k variables for Runge-Kutta calculation of y[n+1]
        elements<T> k1, k2, k3, k4, k5, k6, k7;
// Runge-Kutta algorithm       
//      k1 = h*f(t, y[n])
        k1 = calc_k(stepSize, y[iteration], coeff, accel, *curTime, timeFinal);        
//      k2 = h*f(t+1/5, y[n]+k1*1/5)
        k2 = calc_k(stepSize, y[iteration]+k1*1/5,coeff, accel, *curTime+1/5*stepSize, timeFinal);   
//      k3 = h*f(t+3/10, y[n]+k1*3/40+k2*9/40)
        k3 = calc_k(stepSize, y[iteration]+k1*3/40+k2*9/40,coeff, accel, *curTime+3/10*stepSize, timeFinal);   
//      k4 = h*f(t+4/5, y[n]+k1*44/45+k2*-56/15+k3*32/9)
        k4 = calc_k(stepSize,y[iteration]+k1*44/45+k2*-56/15+k3*32/9,coeff, accel, *curTime+4/5*stepSize, timeFinal);    
//      k5 = h*f(t+8/9, y[n]+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729)
        k5 = calc_k(stepSize, y[iteration]+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729,coeff, accel, *curTime+8/9*stepSize, timeFinal);        
//      k6 = h*f(t, y[n]+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656)
        k6 = calc_k(stepSize, y[iteration]+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656,coeff, accel, *curTime+stepSize, timeFinal);        
//      k7 = h*f(t, y[n]+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84)
        k7 = calc_k(stepSize,y[iteration]+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84,coeff, accel, *curTime+stepSize, timeFinal);  

//      Previous value 
//      v = y[n] + 5179/57600*k1 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7
        *v = y[iteration] + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;     

//      Current value
//      u = y[n] + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
        *u = y[iteration] + k1*(35./384) + k3*(500./1113) + k4*125./192 - k5*2187/6784 + k6*11/84;
}

template <class T> T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, T & stepSize){
//TODO: SC: normTotError might mean Normilize Total Error, but it should be documented. Also what does scale refer to? Why do we need to variables?
    T normTotError, scale;

//      relative error (unitless) 
        elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
        difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);


    //TODO: changing to static cast alters the results slightly

//  square root of the sum of squares of the errors to determine the scale for the time step of the next iteration
//  normalized total error
    normTotError = pow(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2),static_cast<T>(1/2));
    scale = pow((absTol/normTotError),static_cast<T>(1/5));

        return scale;   
}