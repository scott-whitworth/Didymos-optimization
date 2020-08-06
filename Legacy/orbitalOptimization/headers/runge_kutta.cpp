// Didymos-Optimization_Project:
// Last Editor: Mateo and Lauren
// Tasks Completed: 
    // Re-defined the variables going in and out of rkCalc() as well as all functions which call rkCalc() to make intuitive sense.

#include "runge_kutta.h"
#include "acceleration.h" //used for calc_accel() and calc_coast()
#include <iostream> // used for cout
#include <cmath> // used for sine, cosine, and pow functions

template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y_new, 
const T & absTol, coefficients<T> coeff, T & accel, T *gamma,  T *tau, double & lastStep, T *accel_output, T *fuelSpent, const T & wetMass)
{
    // k variables for Runge-Kutta calculation of y[n+1]
    elements<T> k1, k2, k3, k4, k5, k6, k7;

    T curTime = timeInitial; // setting time equal to the start time
    int n=0; // setting the initial iteration number equal to 0
    int minStep=0;

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
    fuelSpent[0]=massFuelSpent;
    elements<T> u, error;
    // Declare the variable for the status of coasting
    bool coast;

    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining deltaT for calc_accel as the stepsize
        T deltaT = stepSize;

        u = y_new[n];

        // defining coast using calc_coast()
        coast = calc_coast(coeff, curTime, timeFinal);
        
        // defining acceleration using calc_accel()
        accel = calc_accel(y_new[n].r,y_new[n].z, NEXT, massFuelSpent, deltaT, coast, wetMass);
        // Record the updated massFuelSpent to the output array
        fuelSpent[n]=massFuelSpent;
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
        if (stepSize>(timeFinal-timeInitial)/100)
        {
            stepSize = (timeFinal-timeInitial)/100;
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
    //std::cout<<"Number of steps: "<<n<<"\n"<<"Min steps :"<<minStep<<"\n";
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
    // Declare the variable for the status of coasting
    bool coast;
    // Declare element for the error of the RK step
    elements<T> error;
    while(curTime<timeFinal) // iterate until time is equal to the stop time
    {
        // defining coast using calc_coast()
        coast = calc_coast(coeff, curTime, timeFinal);

        // defining acceleration using calc_accel()
        accel = calc_accel(y_new.r,y_new.z, NEXT, massFuelSpent, stepSize, coast, static_cast<double>(wetMass));

        //calculate k values
        rkCalc(curTime, timeFinal, stepSize, y_new, coeff, accel, error, k1, k2, k3, k4, k5, k6, k7); 

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        stepSize *= calc_scalingFactor(y_new-error,error,absTol,stepSize)/2;

        // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (stepSize>(timeFinal-timeInitial)/100)
            stepSize = (timeFinal-timeInitial)/100;
        else if (stepSize<((timeFinal-timeInitial)/10000))
            stepSize = (timeFinal-timeInitial)/10000;
        // shorten the last step to end exactly at time final
        if((curTime+stepSize)>timeFinal)
            stepSize = (timeFinal-curTime);

        // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft increases to 1000, so that path is not used for optimization.
        if (sqrt(pow(y_new.r,2)+pow(y_new.z,2))<0.5)
        {
            y_new.r = 1000;
            return;
        }
    }//end of while 
}

template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
T stepSize, elements<T> & y_new, const T & absTol)
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
        rkCalc_Earth(curTime, timeFinal, stepSize, y_new, error, k1, k2, k3, k4, k5, k6, k7);

        //array of time output as t         
        curTime += stepSize;

        //Alter the step size for the next iteration
        //Expected to be negative
        stepSize *= calc_scalingFactor(y_new-error, error,absTol,stepSize)/static_cast <double> (2);

        // The absolute value of step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
        if (-stepSize>(timeFinal-timeInitial)/static_cast <double> (100))
            stepSize = -(timeFinal-timeInitial)/static_cast <double> (100);
        else if (-stepSize<((timeFinal-timeInitial)/static_cast <double> (10000)))
            stepSize = -(timeFinal-timeInitial)/static_cast <double> (10000);

        // shorten the last step to end exactly at time final
        if((curTime+stepSize)<timeInitial)
            stepSize = -(curTime-timeInitial);
    } //end of while 
}

template <class T> void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, coefficients<T> & coeff, const T & accel, 
elements<T> & error, elements<T> k1, elements<T> k2, elements<T> k3, elements<T> k4, elements<T> k5, elements<T> k6, elements<T> k7){
   

    k1 = calc_k(stepSize, y_new, coeff, accel, curTime, timeFinal);      
    k2 = calc_k(stepSize, y_new+k1*static_cast <double> (1)/static_cast <double> (5),coeff, accel, curTime+static_cast <double> (1)/static_cast <double> (5)*stepSize, timeFinal); 
    k3 = calc_k(stepSize, y_new+k1*static_cast <double> (3)/static_cast <double> (40)+k2*static_cast <double> (9)/static_cast <double> (40),coeff, accel, curTime+static_cast <double> (3)/static_cast <double> (10)*stepSize, timeFinal);   
    k4 = calc_k(stepSize, y_new+k1*static_cast <double> (44)/static_cast <double> (45)+k2*static_cast <double> (-56)/static_cast <double> (15)+k3*static_cast <double> (32)/static_cast <double> (9),coeff, accel,curTime+static_cast <double> (4)/static_cast <double> (5)*stepSize, timeFinal); 
    k5 = calc_k(stepSize, y_new+k1*static_cast <double> (19372)/static_cast <double> (6561)+k2*static_cast <double> (-25360)/static_cast <double> (2187)+k3*static_cast <double> (64448)/static_cast <double> (6561)+k4*static_cast <double> (-212)/static_cast <double> (729),coeff, accel, curTime+static_cast <double> (8)/static_cast <double> (9)*stepSize, timeFinal); 
    k6 = calc_k(stepSize, y_new+k1*static_cast <double> (9017)/static_cast <double> (3168)+k2*static_cast <double> (-355)/static_cast <double> (33)+k3*static_cast <double> (46732)/static_cast <double> (5247)+k4*static_cast <double> (49)/static_cast <double> (176)+k5*static_cast <double> (-5103)/static_cast <double> (18656),coeff, accel, curTime+stepSize, timeFinal);  
    k7 = calc_k(stepSize, y_new+k1*static_cast <double> (35)/static_cast <double> (384)+k3*static_cast <double> (500)/static_cast <double> (1113)+k4*static_cast <double> (125)/static_cast <double> (192)+k5*static_cast <double> (-2187)/static_cast <double> (6784)+k6*static_cast <double> (11)/static_cast <double> (84),coeff, accel, curTime+stepSize, timeFinal);  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + k1*static_cast <double> (35)/static_cast <double> (384) + k3*static_cast <double> (500)/static_cast <double> (1113) + k4*static_cast <double> (125)/static_cast <double> (192) - k5*static_cast <double> (2187)/static_cast <double> (6784) + k6*static_cast <double> (11)/static_cast <double> (84);  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    error = k1*static_cast <double> (71)/static_cast <double> (57600) + k3*static_cast <double> (-71)/static_cast <double> (16695) + k4*static_cast <double> (71)/static_cast <double> (1920) - k5*static_cast <double> (17253)/static_cast <double> (339200) + k6*static_cast <double> (22)/static_cast <double> (525) + k7*static_cast <double> (-1)/static_cast <double> (40);    
}

template <class T> void rkCalc_Earth(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error,elements<T> & k1,
elements<T> & k2,elements<T> & k3,elements<T> & k4,elements<T> & k5,elements<T> & k6,elements<T> & k7){
    // Runge-Kutta algorithm      
    //elements<T> y_prev;

    //calc_k multiplies all values by the stepSize internally.
    k1 = calc_k_earth(stepSize, y_new, curTime, timeFinal);        
    k2 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (1)/static_cast <double> (5), curTime+static_cast <double> (1)/static_cast <double> (5)*stepSize, timeFinal);   
    k3 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (3)/static_cast <double> (40)+k2*static_cast <double> (9)/static_cast <double> (40), curTime+static_cast <double> (3)/static_cast <double> (10)*stepSize, timeFinal);   
    k4 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (44)/static_cast <double> (45)+k2*static_cast <double> (-56)/static_cast <double> (15)+k3*static_cast <double> (32)/static_cast <double> (9), curTime+static_cast <double> (4)/static_cast <double> (5)*stepSize, timeFinal);    
    k5 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (19372)/static_cast <double> (6561)+k2*static_cast <double> (-25360)/static_cast <double> (2187)+k3*static_cast <double> (64448)/static_cast <double> (6561)+k4*static_cast <double> (-212)/static_cast <double> (729), curTime+static_cast <double> (8)/static_cast <double> (9)*stepSize, timeFinal);        
    k6 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (9017)/static_cast <double> (3168)+k2*static_cast <double> (-355)/static_cast <double> (33)+k3*static_cast <double> (46732)/static_cast <double> (5247)+k4*static_cast <double> (49)/static_cast <double> (176)+k5*static_cast <double> (-5103)/static_cast <double> (18656), curTime+stepSize, timeFinal);        
    k7 = calc_k_earth(stepSize, y_new+k1*static_cast <double> (35)/static_cast <double> (384)+k3*static_cast <double> (500)/static_cast <double> (1113)+k4*static_cast <double> (125)/static_cast <double> (192)+k5*static_cast <double> (-2187)/static_cast <double> (6784)+k6*static_cast <double> (11)/static_cast <double> (84), curTime+stepSize, timeFinal);  

    //Error 
    //See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    //v = y_new + k1*5179/57600 + k3*7571/16695 + k4*393/640 - k5*92097/339200 + k6*187/2100 + k7*1/40;  

    //New value
    //u = y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    y_new = y_new + k1*static_cast <double> (35)/static_cast <double> (384) + k3*static_cast <double> (500)/static_cast <double> (1113) + k4*static_cast <double> (125)/static_cast <double> (192) - k5*static_cast <double> (2187)/static_cast <double> (6784) + k6*static_cast <double> (11)/static_cast <double> (84);  

    // Error 
    // See the original algorithm by J.R. Dormand and P.J. Prince, JCAM 1980 and its implementation in MATLAB's ode45
    // Dormand-Prince : no error between GPU and CPU
    //y_prev = k1*5179./57600 + k3*7571./16695 + k4*393./640 - k5*92097./339200 + k6*187./2100 + k7*1./40;  
    //error = y_new-y_prev;

    error = k1*static_cast <double> (71)/static_cast <double> (57600) + k3*static_cast <double> (-71)/static_cast <double> (16695) + k4*static_cast <double> (71)/static_cast <double> (1920) - k5*static_cast <double> (17253)/static_cast <double> (339200) + k6*static_cast <double> (22)/static_cast <double> (525) + k7*static_cast <double> (-1)/static_cast <double> (40);    
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
    normTotError = pow(pow(pmError.r,2) + pow(pmError.theta,2) + pow(pmError.z,2) + pow(pmError.vr,2) + pow(pmError.vtheta,2) + pow(pmError.vz,2),(T)0.5);
    scale = pow((absTol/normTotError),(T)1/5);

    return scale;   
}