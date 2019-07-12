#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream>
#include <iomanip> // setprecision(int) 

int main()
{
    // Answers may not align between launchCon() and earthInitial if the startTime and endTime differ

    double startTime = 31557600; // one year (s)
    double endTime = 94672800; // three years (s)
    double timeRes = 3600; // minute resolution
    double testPoint = startTime+330; // testing (s)

    EarthInfo launchCon(startTime,endTime,timeRes);

    elements<double> result;
    result = launchCon.getCondition(testPoint);

    /*for(int i = 0 ; i < launchCon.getTolData(); i++){
        elements<double> result;
        result = launchCon.getCondition(startTime + (i * timeRes) + 25);

        std::cout << i << "th possible output: " << result << std::endl;
    }*/

}