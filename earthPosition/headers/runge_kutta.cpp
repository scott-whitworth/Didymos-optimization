#include "runge_kutta.h"
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions


template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> & y, const T & absTol)
{
    // Set the first element of the solution vector to the conditions of earth on impact date (Oct. 5, 2022)
    y = y0;
    // k variables for Runge-Kutta calculation of y for earth's initial position (launch date)
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    T curTime = timeFinal; // setting time equal to the start time

    while(curTime>timeInitial) // iterates in reverse
    {
        elements<T> v;

        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y, v, y);

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        //Expected to be negative
        stepSize *= calc_scalingFactor(v,y-v,absTol,stepSize)/50;

        // The absolute value of step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (-stepSize>(timeFinal-timeInitial)/100)
            stepSize = -(timeFinal-timeInitial)/100;
        else if (-stepSize<((timeFinal-timeInitial)/100000))
            stepSize = -(timeFinal-timeInitial)/100000;
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)<timeInitial)
            stepSize = -(curTime-timeInitial);
    }//end of while 
}

template <class T> void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> y, elements<T> & v, elements<T> & u){
    // Runge-Kutta algorithm      
    elements<T> k1, k2, k3, k4, k5, k6, k7; 
    //k1 = h*f(t, y)
    k1 = calc_k(stepSize, y, curTime, timeFinal);        
    //k2 = h*f(t+1/5, y+k1*1/5)
    k2 = calc_k(stepSize, y+k1*1/5, curTime+1/5*stepSize, timeFinal);   
    //k3 = h*f(t+3/10, y+k1*3/40+k2*9/40)
    k3 = calc_k(stepSize, y+k1*3/40+k2*9/40, curTime+3/10*stepSize, timeFinal);   
    //k4 = h*f(t+4/5, y+k1*44/45+k2*-56/15+k3*32/9)
    k4 = calc_k(stepSize,y+k1*44/45+k2*-56/15+k3*32/9, curTime+4/5*stepSize, timeFinal);    
    //k5 = h*f(t+8/9, y+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729)
    k5 = calc_k(stepSize, y+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729, curTime+8/9*stepSize, timeFinal);        
    //k6 = h*f(t, y+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656)
    k6 = calc_k(stepSize, y+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656, curTime+stepSize, timeFinal);        
    //k7 = h*f(t, y+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84)
    k7 = calc_k(stepSize,y+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84, curTime+stepSize, timeFinal);  

    // Previous value 
    //v = y + 5179/57600*k1 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7
    v = y + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;     

    //Current value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    u = y + k1*(35./384) + k3*(500./1113) + k4*125./192 - k5*2187/6784 + k6*11/84;  
}

template <class T> T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, T & stepSize)
{
    // relative total error is the total error of all coponents of y which is used in scale.
    // scale is used to determine the next step size.
    T normTotError, scale;

    // relative error (unitless) 
    elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
    difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);

    // square root of sum of squares of the error from the 6 elements to determine the scale for the time step of the next iteration
    normTotError = pow(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2),(T)1/2);
    scale = pow((absTol/normTotError),(T)1/5);

    return scale;   
}