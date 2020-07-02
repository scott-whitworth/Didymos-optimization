#ifndef earthInfo_cpp
#define earthInfo_cpp

#include "earthInfo.h"
    
#include <iostream>  // cout
#include <iomanip> // setprecision(int)  


EarthInfo::EarthInfo(const double & beginTime, const double & stopTime, const double & timeAcc, const cudaConstants* cConstants) {
    ////////////////////////////////////
    ///Setting up initial information///
    ////////////////////////////////////

    startTime = beginTime; // Starting time (s)
    endTime = stopTime; // Ending time (s)
    timeRes = timeAcc; // Time resolution (s)
    tolData = ((endTime-startTime)/timeRes) + 1; // Total Number of Data points (in hours) in earthCon, plus one for the last 'section'
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Alocate memory for earthCon, one entry for every data point
    earthCon = new elements<double> [tolData];

    // Assigning the position of the earth at impact to variable earth. Passed into earthInitial_incremental and rk4Reverse.
    elements<double> earth = elements<double>(cConstants->r_fin_earth, cConstants->theta_fin_earth, cConstants->z_fin_earth, cConstants->vr_fin_earth, cConstants->vtheta_fin_earth, cConstants->vz_fin_earth);

    // Get the actual initial position for the time frame chosen.
    earth = earthInitial_incremental(0, startTime, earth, cConstants);

    // Setting the first element of earthCon to be equal to the earth conditions at impact
    //earthCon[0]=earthInitial_incremental(startTime,startTime,earth);
    earthCon[0] = earth;

    // Shows progress of earth position calculations before the optimization in cuda can occur.
    std::cout << std::endl << "calculating earth positions for the trip time range" << std::endl;
    std::cout << "          10 20 30 40 50 60 70 80 90 100" << std::endl;
    std::cout << "progress:[";

    for(int i=1; i < tolData; i++) { 
        // Calculates earth's condition at each point (time) from the previously calculated point.
        earthCon[i] = earthInitial_incremental(calc_time(i)-timeRes, calc_time(i), earth, cConstants); // Obtaining conditions of the earth

        // Filling progress bar
        if ((i % (tolData/30)) == 0) {
            std::cout << ">";
        }

        // Sets earth equal to the conditions calculated for a given time.
        earth = earthCon[i];
    }
    // Closing progress bar
    std::cout << "]\n\n";
}


elements<double> EarthInfo::getCondition(const double & currentTime) {
    // Defining variables //
    elements<double> lower;
    elements<double> upper;
    int index;

    // Setting index equal to the location of data based on time
    index = calcIndex(currentTime); 
    // The lower index weight is equal to the index
    lower = earthCon[index];
    // The upper index weight is equal to the index plus 1
    upper = earthCon[index + 1];
    
    // Weights are used for interpolation
    double lowerWeight = 1 - ((currentTime-calc_time(index))/timeRes);
    double upperWeight = 1 - ((calc_time(index+1)-currentTime)/timeRes);
    elements<double> result = interpolate(lower,upper,lowerWeight,upperWeight);
/*    
    if(result.r<=0.1)
    {
        std::cout<<"NaN incoming \n" << result;
        std::cout<<"LowerW : "<<lowerWeight<<"  UpperW : "<<upperWeight << "\n";
        std::cout<<"Timne 1 : "<<calc_time(index)<<"  Time 2 : "<<calc_time(index+1) << "\n";
    }
*/
    return result;
}

int EarthInfo::calcIndex(const double & currentTime) {
   return static_cast<int>((currentTime-startTime)/timeRes);

}

double EarthInfo::calc_time(const int & currentIndex) {
    return startTime + (currentIndex*timeRes);
}

int EarthInfo::getTolData() {
    return tolData;
}

elements<double> EarthInfo::interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight) {
    return (lower*lowerWeight)+(upper*upperWeight);
}

EarthInfo::~EarthInfo() {
    delete [] earthCon;
}

elements<double> earthInitial_incremental(double timeInitial, double tripTime, const elements<double> & earth, const cudaConstants * cConstants) {
  // Time step
  double deltaT; 

  // Initial guess for time step, cannot be greater than the time resolution.
  deltaT = -(tripTime - timeInitial)/static_cast <double> (60); 

  // Declaring the solution vector.
  elements<double> yp;

  // Calculates the earth's launch date conditions based on timeFinal minus the optimized trip time.
  rk4Reverse(timeInitial,tripTime,earth,deltaT,yp, cConstants->rk_tol);
 
  return yp;
}

#endif