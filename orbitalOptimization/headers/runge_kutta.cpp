#include "runge_kutta.h"
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y, 
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
    y[0] = y0;
    times[0]=timeInitial;
    // array of gamma for binary output
    gamma[0] = calc_gamma(coeff,timeInitial, timeFinal);
    // array of tau for binary output
    tau[0] = calc_tau(coeff,timeInitial, timeFinal); 
    // array of acceleration for binary output
    accel_output[0] = calc_accel(y[0].r, NEXT, massFuelSpent, stepSize, calc_coast(coeff, curTime, timeFinal), wetMass);

    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining deltaT for calc_accel as the stepsize
        T deltaT = stepSize;

        // defining coast using calc_coast()
        T coast = calc_coast(coeff, curTime, timeFinal);

        // defining acceleration using calc_accel()
        accel = calc_accel(y[n].r, NEXT, massFuelSpent, deltaT, coast, wetMass);

        // to hold previous and  current  values
        elements<T> v, u;
        
        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y[n], coeff, accel, v, u);

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
        stepSize *= calc_scalingFactor(v,u-v,absTol,stepSize);

        //The step size cannot exceed the total time divided by 10 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/2)
        {
            stepSize = (timeFinal-timeInitial)/2;
            maxStep++;
        }
        else if (stepSize<((timeFinal-timeInitial)/1000))
        {
            stepSize = (timeFinal-timeInitial)/1000;
            minStep++;
        }
        if((curTime+stepSize)>timeFinal)
            stepSize = (timeFinal-curTime);


        //Calculates the y[n] for the next round of calculations
        y[n+1] = u;   
        n++;
    }//end of while 
    std::cout<<massFuelSpent<<std::endl;
    lastStep = n;
    std::cout<<"Number of steps: "<<n<<"\n"<<"Min steps :"<<minStep<<"\n"<<"Max steps: "<<maxStep<<"\n";
}



template <class T> void rk4Simple(const T & timeInitial, const T & timeFinal, const elements<T> & y0,
T stepSize, elements<T> & y, const T & absTol, coefficients<T> coeff, T & accel, const T & wetMass)
{
    // Set the first element of the solution vector to the initial conditions of the spacecraft
    y = y0;
    // k variables for Runge-Kutta calculation of y based off the spacecraft's final state
    elements<T> k1, k2, k3, k4, k5, k6, k7;
    T curTime = timeInitial; // setting time equal to the start time

    // corresponds NEXT thruster to type 1 in thruster.h
    thruster<T> NEXT = thruster<T>(1);

    //mass of fuel expended (kg)
    //set to 0 initially
    T massFuelSpent =0;

    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining deltaT for calc_accel as the stepsize
        T deltaT = stepSize;

        // defining coast using calc_coast()
        T coast = calc_coast(coeff, curTime, timeFinal);

        // defining acceleration using calc_accel()
        accel = calc_accel(y.r, NEXT, massFuelSpent, deltaT, coast, wetMass);
        elements<T> v;


        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y, coeff, accel, v, y); 

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(v,y-v,absTol,stepSize);

        // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/2)
            stepSize = (timeFinal-timeInitial)/2;
        else if (stepSize<((timeFinal-timeInitial)/1000))
            stepSize = (timeFinal-timeInitial)/1000;
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)>timeFinal)
            stepSize = (timeFinal-curTime);

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft increases to 1000, so that path is not used for optimization.
        if (y.r<0.5)
        {
            y.r = 1000;
        }
    }//end of while 
}

template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> & y, const T & absTol, coefficients<T> coeff, const T & accel)
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
        rkCalc(curTime, timeFinal, stepSize, y, coeff, accel, v, y);

        //array of time output as t         
        curTime = curTime + stepSize;

        //Alter the step size for the next iteration
        //Expected to be negative
        stepSize *= calc_scalingFactor(v,y-v,absTol,stepSize);

        // The absolute value of step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (-stepSize>(timeFinal-timeInitial)/2)
            stepSize = -(timeFinal-timeInitial)/2;
        else if (-stepSize<((timeFinal-timeInitial)/1000))
            stepSize = -(timeFinal-timeInitial)/1000;
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)<timeInitial)
            stepSize = -(curTime-timeInitial);
    }//end of while 
}

template <class T> void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> y, coefficients<T> & coeff, const T & accel, elements<T> & v, elements<T> & u){
    // Runge-Kutta algorithm      
    elements<T> k1, k2, k3, k4, k5, k6, k7; 
    //k1 = h*f(t, y)
    k1 = calc_k(stepSize, y, coeff, accel, curTime, timeFinal);        
    //k2 = h*f(t+1/5, y+k1*1/5)
    k2 = calc_k(stepSize, y+k1*1/5,coeff, accel, curTime+1/5*stepSize, timeFinal);   
    //k3 = h*f(t+3/10, y+k1*3/40+k2*9/40)
    k3 = calc_k(stepSize, y+k1*3/40+k2*9/40,coeff, accel, curTime+3/10*stepSize, timeFinal);   
    //k4 = h*f(t+4/5, y+k1*44/45+k2*-56/15+k3*32/9)
    k4 = calc_k(stepSize,y+k1*44/45+k2*-56/15+k3*32/9,coeff, accel, curTime+4/5*stepSize, timeFinal);    
    //k5 = h*f(t+8/9, y+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729)
    k5 = calc_k(stepSize, y+k1*19372/6561+k2*-25360/2187+k3*64448/6561+k4*-212/729,coeff, accel, curTime+8/9*stepSize, timeFinal);        
    //k6 = h*f(t, y+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656)
    k6 = calc_k(stepSize, y+k1*9017/3168+k2*-355/33+k3*46732/5247+k4*49/176+k5*-5103/18656,coeff, accel, curTime+stepSize, timeFinal);        
    //k7 = h*f(t, y+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84)
    k7 = calc_k(stepSize,y+k1*35/384+k3*500/1113+k4*125/192+k5*-2187/6784+k6*11/84,coeff, accel, curTime+stepSize, timeFinal);  

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

    //TODO: changing to static cast alters the results slightly

    return scale;   
}