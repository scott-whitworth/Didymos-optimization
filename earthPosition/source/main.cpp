#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream>
#include <iomanip> // setprecision(int) 

int main()
{
    double startTime = 31557600; // one year (s)
    double endTime = 94672800; // three years (s)
    double timeRes = 3600; // minute resolution
    double testPoint = startTime+330; // testing (s)

    EarthInfo launchCon(startTime,endTime,timeRes);

    elements<double> result;
    result = launchCon.getCondition(testPoint);


    std::ofstream output;
  
    output.open ("getEarthInfo.bin", std::ios::binary); 
    double timetime =0;
    for(int i = 0 ; i < launchCon.getTolData(); i++)
    {
        timetime = startTime + (i * timeRes);
        elements<double> result;
        result = launchCon.getCondition(timetime);
        output.write((char*)&result, sizeof (elements<double>));
        output.write((char*)&timetime, sizeof (double));

    }
    output.close();
}