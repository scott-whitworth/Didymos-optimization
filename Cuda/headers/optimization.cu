//Didymos-Optimization_Project:
//Last Editor: Ben
//Tasks Completed: 
    //Put for loop in main to call new optimize() function


#include "optimization.h" 
#include "constants.h" //used for wetMass
#include "orbitalMotion.h" //used for trajectory() and trajectoryPrint()
#include "earthInfo.h"

#include <iostream> // cout
#include <iomanip> //used for setw(), sets spaces between values output
#include <time.h> //for seeding the random number generator
#include <random>
#include <chrono>


#include "runge_kuttaCUDA.cuh" //for testing rk4simple


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
    //int blockThreadNums[] = {32, 64, 192, 256, 384, 512, 768, 1024};
    int blockThreadNums[] = {32};
    //int threadNums[] = {100, 500, 1000, 2000, 3000, 4000, 5000};
    int threadNums[] = {2000};

    std::ofstream efficiencyGraph;
    efficiencyGraph.open("efficiencyGraph.csv");

    double calcPerS;

    for(int i = 0; i < std::size(blockThreadNums); i++){
        for(int j = 0; j < std::size(threadNums); j++){    
            blockThreads = blockThreadNums[i];
            numThreads = threadNums[j];
            std::cout << std::endl << "testing optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads" << std::endl;
            //calcPerS = optimize(numThreads, blockThreads);
            optimize(numThreads, blockThreads);
            //efficiencyGraph << blockThreads << "," << numThreads << "," << calcPerS  << "\n";
        }
    }

    efficiencyGraph.close();
    
    delete launchCon;
    return 0;
}