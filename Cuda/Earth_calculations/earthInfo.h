#ifndef earthInfo_h
#define earthInfo_h

#include "../Motion_Eqns/elements.h"
#include "../Config_Constants/config.h"

class EarthInfo {
    private:

    // Holds all of the earth conditions for a given time range
    elements<double> *earthCon;

    // First time point for a given time span
    double startTime;

    // Last time point for a given time span
    double endTime;

    // The resolution of data points (ex: 3600 = hours, 60 = minutes, 1 = seconds)
    double timeRes;

    // The total amount of data points for a run. Calculated by time span divided by timeRes.
    int tolData;

    // Takes in a time and outputs a corresponding index (location of data).
    int calcIndex(const double & currentTime);

    // Takes in an index and outputs the time corresponded to that index.
    double calc_time(const int & currentIndex);

    // Takes in the lower index, upper index, and their weights in order to calculate the Earth's position for a time between two index.
    elements<double> interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:

    // Used to initialize the earth calculation data. Fills up earthCon.
    EarthInfo(const double & beginTime, const double & stopTime, const double & timeAcc, cudaConstants& cConstants);

    // Returns the interpolated conditions of earth for a given time input.
    elements<double> getCondition(const double & currentTime);

    // Returns the total amount of data for a run with a given time span and resolution.
    int getTolData();

    // Clears dynamic memory after each run
    ~EarthInfo();

};

#include "earthInfo.cpp"

// Makes a global variable for launchCon (called in optimization.cu)
EarthInfo *launchCon;

#endif 