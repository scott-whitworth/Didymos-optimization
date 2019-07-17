#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream> // used for binary output
#include <iomanip> // setprecision(int) 

int main()
{
    // Define variables to be passed into EarthInfo
    //double startTime = 47304000; // 1.5 year (s)
    //double endTime = 78840000; // 2.5 years (s)
    double startTime = 0; // 0 year (s)
    double endTime = 94608000; // 3 years (s)
    double timeRes = 3600; // seconds to hours conversion

    // initializes EarthInfo
    EarthInfo launchCon(startTime,endTime,timeRes);


    // output to binary file
    std::ofstream output;
    output.open ("getEarthInfo.bin", std::ios::binary); 
    double timetime =0.;
    for(int i = 0 ; i <= launchCon.getTolData(); i++)
    {
        timetime = startTime + (i * timeRes);//creating a timeline in seconds, with a tick every hour
        elements<double> result;
        result = launchCon.getCondition(timetime);//updating every 1 hour
        output.write((char*)&result, sizeof (elements<double>));
        output.write((char*)&timetime, sizeof (double));
    }
    output.close();
}