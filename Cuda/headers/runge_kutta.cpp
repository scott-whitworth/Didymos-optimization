//Didymos-Optimization_Project:
//Last Editor: Ben and Lauren
//Tasks Completed: 
    //Functionalized rkCaLc() which is called by all three of the runge-kutta functions.
    //Added the z component to the calcAccel() function calls

#include "runge_kutta.h"
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y_new, 
const T & absTol, coefficients<T> coeff, T & accel, T *gamma,  T *tau, int & lastStep, T *accel_output, const T & wetMass)
{
    // k variables for Runge-Kutta calculation of y[n+1]
    elements<T> k1, k2, k3, k4, k5, k6, k7;

    T curTime = timeInitial; // setting time equal to the start time
    int n=0; // setting the initial iteration number equal to 0
    int minStep=0;
    int maxStep=0;

    // corresponds NEXT thruster to type 1 in thruster.h
    thruster<T> NEXT = thruster<T>(1);

    //mass of fuel expended (kg)
    //set to 0 initially
    T massFuelSpent =0;

    // Set the first element of the solution vector to the initial conditions
    y_new[0] = y0;
    times[0]=timeInitial;
    // array of gamma for binary output
    gamma[0] = calc_gamma(coeff,timeInitial, timeFinal);
    // array of tau for binary output
    tau[0] = calc_tau(coeff,timeInitial, timeFinal); 
    // array of acceleration for binary output
    accel_output[0] = calc_accel(y_new[0].r,y_new[0].z, NEXT, massFuelSpent, stepSize, calc_coast(coeff, curTime, timeFinal), wetMass);

    elements<T> u, error;

    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining deltaT for calc_accel as the stepsize
        T deltaT = stepSize;

        u = y_new[n];

        // defining coast using calc_coast()
        bool coast = calc_coast(coeff, curTime, timeFinal);
        
        // defining acceleration using calc_accel()
        accel = calc_accel(y_new[n].r,y_new[n].z, NEXT, massFuelSpent, deltaT, coast, wetMass);
        
        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, u, coeff, accel, error, k1, k2, k3, k4, k5, k6, k7);

        //array of time output as t         
        curTime += stepSize;
        //Time of iteration is set to the previous time plus the step size used within that iteration
        times[n+1]=curTime;
        //array of gamma for binary output
        gamma[n+1] = calc_gamma(coeff,curTime, timeFinal);
        //array of tau for binary output
        tau[n+1] = calc_tau(coeff,curTime, timeFinal);  
        //array of accel for binary output
        accel_output[n+1]= accel;


        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(u-error,error,absTol,stepSize)/2;

        //The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/1000)
        {
            stepSize = (timeFinal-timeInitial)/1000;
            maxStep++;
        }
        else if (stepSize<((timeFinal-timeInitial)/10000))
        {
            stepSize = (timeFinal-timeInitial)/10000;
            minStep++;
        }
        
        if((curTime+stepSize)>timeFinal)
            stepSize = (timeFinal-curTime);


        //Calculates the y[n] for the next round of calculations
        y_new[n+1] = u;   
        n++;
    }//end of while 
    lastStep = n;
    std::cout<<"Number of steps: "<<n<<"\n"<<"Min steps :"<<minStep<<"\n"<<"Max steps: "<<maxStep<<"\n";
}



template <class T> void rk4Simple(const T & timeInitial, const T & timeFinal, const elements<T> & y0,
T stepSize, elements<T> & y_new, const T & absTol, coefficients<T> coeff, T & accel, const T & wetMass)
{
    // Set the first element of the solution vector to the initial conditions of the spacecraft
    y_new = y0;
    // k variables for Runge-Kutta calculation of y based off the spacecraft's final state
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    T curTime = timeInitial; // setting time equal to the start time

    // corresponds NEXT thruster to type 1 in thruster.h
    thruster<T> NEXT = thruster<T>(1);

    //mass of fuel expended (kg)
    //set to 0 initially
    T massFuelSpent =0;

    elements<T> error;

    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining coast using calc_coast()
        bool coast = calc_coast(coeff, curTime, timeFinal);

        // defining acceleration using calc_accel()
        accel = calc_accel(y_new.r,y_new.z, NEXT, massFuelSpent, stepSize, coast, wetMass);
        //accel = 0.;

        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y_new, coeff, accel, error, k1, k2, k3, k4, k5, k6, k7); 

        //array of time output as t         
        curTime += stepSize;

        
        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol,stepSize)/2;

        // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/1000)
            stepSize = (timeFinal-timeInitial)/1000;
        else if (stepSize<((timeFinal-timeInitial)/10000))
            stepSize = (timeFinal-timeInitial)/10000;
        
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)>timeFinal)
            stepSize = (timeFinal-curTime);
        

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft increases to 1000, so that path is not used for optimization.
        if (y_new.r<0.5)
        {
            y_new.r = 1000;
        }
    }//end of while 
}

template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> & y_new, const T & absTol, coefficients<T> coeff, const T & accel)
{
    // Set the first element of the solution vector to the conditions of earth on impact date (Oct. 5, 2022)
    y_new = y0;
    // k variables for Runge-Kutta calculation of y for earth's initial position (launch date)
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    elements<T> error;
    T curTime = timeFinal; // setting time equal to the start time

    while(curTime>timeInitial) // iterates in reverse
    {

        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y_new, coeff, accel, error, k1, k2, k3, k4, k5, k6, k7);

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        //Expected to be negative
        stepSize *= calc_scalingFactor(y_new-error, error,absTol,stepSize)/2;

        // The absolute value of step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (-stepSize>(timeFinal-timeInitial)/1000)
            stepSize = -(timeFinal-timeInitial)/1000;
        else if (-stepSize<((timeFinal-timeInitial)/10000))
            stepSize = -(timeFinal-timeInitial)/10000;
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)<timeInitial)
            stepSize = -(curTime-timeInitial);
    }//end of while 
}

template <class T> __host__ __device__ void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, coefficients<T> & coeff, const T & accel, 
elements<T> & error, elements<T> k1, elements<T> k2, elements<T> k3, elements<T> k4, elements<T> k5, elements<T> k6, elements<T> k7){
   
    //coefficients from MATLAB's implementation of ode45
    //Our calculation of k has the time step built into it (see motion_equations.cpp)
    k1 = calc_k(stepSize, y_new,                                                                                coeff, accel, curTime,              timeFinal);      
    k2 = calc_k(stepSize, y_new+k1*1/5,                                                                         coeff, accel, curTime+1/5*stepSize, timeFinal); 
    k3 = calc_k(stepSize, y_new+k1*3/40+k2*9/40,                                                                coeff, accel, curTime+3/10*stepSize,timeFinal);   
    k4 = calc_k(stepSize, y_new+k1*44/45+k2*-56/15+k3*32/9,                                                     coeff, accel, curTime+4/5*stepSize, timeFinal); 
    k5 = calc_k(stepSize, y_new+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729,                         coeff, accel, curTime+8/9*stepSize, timeFinal); 
    k6 = calc_k(stepSize, y_new+k1*9017/3168 +k2*-355/33    +k3*46732/5247+k4*49/176 +k5*-5103/18656,           coeff, accel, curTime+stepSize,     timeFinal);  
    k7 = calc_k(stepSize, y_new+k1*35/384                   +k3*500/1113  +k4*125/192+k5*-2187/6784+k6*11/84,   coeff, accel, curTime+stepSize,     timeFinal);  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + k1*(35./384) + k3*(500./1113) + k4*125./192 - k5*2187./6784 + k6*11./84;  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    error = k1*71./57600 + k3*-71./16695 + k4*71./1920 - k5*17253./339200 + k6*22./525 + k7*-1./40;  
    //error = k1*71./57600 + k3*-71./16695 + k4*71./1920 - k5*17253./339200 + k6*22./525;  
    
}

template <class T> __host__ __device__ T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, T & stepSize)
{
    // relative total error is the total error of all coponents of y which is used in scale.
    // scale is used to determine the next step size.
    T normTotError, scale;

    // relative error (unitless) 
    elements<T> pmError(difference.r/previous.r, difference.theta/previous.theta, difference.z/previous.z, 
    difference.vr/previous.vr,  difference.vtheta/previous.vtheta, difference.vz/previous.vz);

    //elements<T> pmError(previous.r, previous.theta, previous.z, previous.vr,  previous.vtheta, previous.vz);

    // square root of sum of squares of the error from the 6 elements to determine the scale for the time step of the next iteration
    normTotError = pow(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2),(T)1/2);
    scale = pow((absTol/normTotError),(T)1/5);

    //TODO: changing to static cast alters the results slightly

    return scale;   
}