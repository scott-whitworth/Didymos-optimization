//Didymos-Optimization_Project:
//Last Editor: Ben
//Tasks Completed: 
    //Put for loop in main to call new optimize() function

#include "../constants.h" //used for wetMass
#include "../Earth_calculations/orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include "../Earth_calculations/earthInfo.h"
#include "../Runge_Kutta/runge_kuttaCUDA.cuh" //for testing rk4simple

#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values output
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>


int main ()
{
    //////////////////////////////////////////////////////////////////////////////////
    //Global variable needs to be initialized

    // Define variables to be passed into EarthInfo
    double startTime = 31536000; // 1.0 year (s)
    double endTime = 63072000; // 2.0 years (s)
    double timeRes = 3600; // position of earth is calculated for every minute

    // initializes EarthInfo
    launchCon = new EarthInfo(startTime, endTime, timeRes);
    ////////////////////////////////////////////////////////////////////////////////////

    int blockThreads = 0;
    int numThreads = 0;
    int blockThreadNums[] = {32};
    int threadNums[] = {300};

    std::ofstream efficiencyGraph;
    efficiencyGraph.open("efficiencyGraph.csv");

    double calcPerS;

    for(int i = 0; i < std::size(blockThreadNums); i++){
        for(int j = 0; j < std::size(threadNums); j++){    
            blockThreads = blockThreadNums[i];
            numThreads = threadNums[j];
            std::cout << std::endl << "testing optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads" << std::endl;
            calcPerS = optimize(numThreads, blockThreads);
            efficiencyGraph << blockThreads << "," << numThreads << "," << calcPerS  << "\n";
        }
    }

    efficiencyGraph.close();
    
    delete launchCon;
    return 0;
}