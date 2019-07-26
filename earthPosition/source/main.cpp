#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream> // used for binary output
#include <iomanip> // setprecision(int) 
#include <chrono>

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // Define variables to be passed into EarthInfo
    double startTime = 31536000; // 1.0 year (s)
    double endTime = 63072000; // 2.0 years (s)
    double timeRes = 3600; // position of earth is calculated for every minute

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

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

}