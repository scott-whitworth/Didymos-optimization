#ifndef earthInfo_h
#define earthInfo_h

#include "../Motion_Eqns/elements.h"

class EarthInfo
{

    private:
    elements<double> *earthCon;
    double startTime;
    double endTime;
    double timeRes;
    int tolData;

    int calcIndex(const double & currentTime);
    double calc_time(const int & currentIndex);
    elements<double> interpolate(const elements<double> & lower,const elements<double> & upper,const double & lowerWeight,const double & upperWeight);

    public:
    EarthInfo(const double & beginTime, const double & stopTime, const double & timeAcc);
    elements<double> getCondition(const double & currentTime);
    int getTolData();
    ~EarthInfo();

};

#include "earthInfo.cpp"

EarthInfo *launchCon;

#endif 