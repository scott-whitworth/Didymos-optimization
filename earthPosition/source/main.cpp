#include "orbitalMotion.h"
#include "elements.h"
#include "earthInfo.h"

#include <iostream>

int main()
{
    double startTime = 31557600; // one year (s)
    double endTime = 94672800; // three years (s)
    double timeRes = 3600; // minute resolution
    double testPoint = startTime+150000; // testing (s)

    EarthInfo launchCon(startTime,endTime,timeRes);

    launchCon.getCondition(testPoint);

    launchCon.~EarthInfo();
}