#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream>
#include <iomanip> // setprecision(int) 

int main()
{
    double startTime = 31557600; // one year (s)
    double endTime = 94672800; // three years (s)
    double timeRes = 3600; // seconds to hours conversion

    EarthInfo launchCon(startTime,endTime,timeRes);

    std::ofstream output;
  
    output.open ("getEarthInfo.bin", std::ios::binary); 
    double timetime =0., timeJulianDays=0.;
    for(int i = 0 ; i < launchCon.getTolData(); i++)
    {
        timetime = startTime + (i * timeRes);//creating a timeline in seconds, with a tick every hour
        elements<double> result;
        result = launchCon.getCondition(timetime);//updating every 1 hour
        output.write((char*)&result, sizeof (elements<double>));
        timeJulianDays = 2459857.500000000-(timetime/(24*3600));//creating a new output time in Julian days to match JPL data JULIAN TIME
        output.write((char*)&timetime, sizeof (double));
    }
    output.close();
}