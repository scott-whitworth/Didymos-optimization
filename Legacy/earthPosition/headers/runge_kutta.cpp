#include "runge_kutta.h"
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions


template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> & y_new, const T & absTol)
{

    //Called in orbitalMotion.cpp as 'rk4Reverse(timeInitial,tripTime,earth,deltaT,yp,absTol);'
    //--> timeInitial=0.
    //--> timeFinal=tripTime.
    //--> y0 contains the conditions of the earth on the asteroid impact date.
    //--> stepSize is <0 and contains the initial guess for deltaT based on MAX_NUMSTEPS.
    //--> yp is merely a return vector.

    // Set the first element of the solution vector to the conditions of earth on impact date (Oct. 5, 2022)
    y_new = y0;
    // k variables for Runge-Kutta calculation of y for earth's initial position (launch date)
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    elements<T> error; 
    T curTime = timeFinal; // setting time equal to the start time
    int count;

    while(curTime>timeInitial) // iterates in reverse
    {
        count++;
        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y_new, error,k1, k2, k3, k4, k5, k6, k7);

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size (<0) for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol,stepSize)/2;

        //Set limits on the stepSize returned from the previous step
        if (-stepSize>(timeFinal-timeInitial)/100)
            stepSize = -(timeFinal-timeInitial)/100;//Maximum allowed (absolute) value
        else if (-stepSize<((timeFinal-timeInitial)/1000))
            stepSize = -(timeFinal-timeInitial)/1000;//Minimum allowed (absolute) value
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)<timeInitial)
            stepSize = -(curTime-timeInitial);
    }//end of while 
    std::cout << "number of iterations: " << count << std::endl;
}

template <class T> void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error,elements<T> & k1,
elements<T> & k2,elements<T> & k3,elements<T> & k4,elements<T> & k5,elements<T> & k6,elements<T> & k7){
    // Runge-Kutta algorithm      
    elements<T> y_prev;

    //calc_k multiplies all values by the stepSize internally.
    k1 = calc_k(stepSize, y_new, curTime, timeFinal);        
    k2 = calc_k(stepSize, y_new+k1*1/5, curTime+1/5*stepSize, timeFinal);   
    k3 = calc_k(stepSize, y_new+k1*3/40+k2*9/40, curTime+3/10*stepSize, timeFinal);   
    k4 = calc_k(stepSize, y_new+k1*44/45+k2*-56/15+k3*32/9, curTime+4/5*stepSize, timeFinal);    
    k5 = calc_k(stepSize, y_new+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729, curTime+8/9*stepSize, timeFinal);        
    k6 = calc_k(stepSize, y_new+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656, curTime+stepSize, timeFinal);        
    k7 = calc_k(stepSize, y_new+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84, curTime+stepSize, timeFinal);  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    //v = y_new + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + k1*(35./384) + k3*(500./1113) + k4*125./192 - k5*2187./6784 + k6*11./84;  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    error = y_new-y_prev;
}

template <class T> T calc_scalingFactor(const elements<T> & previous, const elements<T> & error, const T & absTol, T & stepSize)
{
    // relative total error is the total error of all coponents of y which is used in scale.
    // scale is used to determine the next step size.
    T normTotError, scale;

    // relative error (unitless) 
    elements<T> pmError(error.r/previous.r, error.theta/previous.theta, error.z/previous.z, 
    error.vr/previous.vr,  error.vtheta/previous.vtheta, error.vz/previous.vz);

    // square root of sum of squares of the error from the 6 elements to determine the scale for the time step of the next iteration
    normTotError = pow(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2),(T)1/2);
    scale = pow((absTol/normTotError),(T)1/5);

    return scale;   
}