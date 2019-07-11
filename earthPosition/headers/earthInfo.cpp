#include "earthInfo.h"
#include "orbitalMotion.h"
    


EarthInfo::EarthInfo(const double & beginTime, const double & stopTime, const double & timeAcc)
{
    //elements<double> *earthCon;
    startTime = beginTime;
    endTime = stopTime;
    timeRes = timeAcc;
    tolData = (endTime-startTime)/timeRes;

    earthCon = new elements<double> [tolData]();
    
    //TODO: iterative runge kutta to speed up execution speed, reuse previous solution.

    for(int i=0; i<tolData; i++)
    {
        earthCon[i]=earthInitial(calc_time(i));
    }
}

elements<double> EarthInfo::interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight)
{
    return (lower*lowerWeight)+(upper*upperWeight);
}

elements<double> EarthInfo::getCondition(const double & currentTime)
{
    elements<double> lower;
    elements<double> upper;
    int index;
    index = calcIndex(currentTime);
    lower = earthCon[index];
    upper = earthCon[index + 1];
    double lowerWeight = (currentTime-calc_time(index))/timeRes;
    double upperWeight = (calc_time(index+1)-currentTime)/timeRes;

    return interpolate(lower,upper,lowerWeight,upperWeight);
}

int EarthInfo::calcIndex(const double & currentTime)
{
   return static_cast<int>((currentTime-startTime)/timeRes);

}

double EarthInfo::calc_time(const int & currentIndex)
{
    return startTime + (currentIndex*timeRes);
}

EarthInfo::~EarthInfo()
{
    delete [] earthCon;
}