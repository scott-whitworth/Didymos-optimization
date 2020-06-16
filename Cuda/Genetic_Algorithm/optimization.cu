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


int main () {
    // display GPU properties and ensure we are using the right one
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Device Number: 0 \n";
    std::cout << "- Device name: " << prop.name << std::endl;
    cudaSetDevice(0);
    
    // pre-calculate a table of Earth's position within possible mission time range
    //----------------------------------------------------------------
    // Define variables to be passed into EarthInfo
    double startTime = 15778800; // 0.5 year (s)
    double endTime = 78894000; // 2.5 years (s)
    double timeRes = 3600; // (s) position of earth is calculated for every hour

    launchCon = new EarthInfo(startTime, endTime, timeRes); // a global variable to hold Earth's position over time
    //----------------------------------------------------------------

    int blockThreads = 32;
    int numThreads = 2880; // the number of cores on a Tesla k40
    //int numThreads = 1920; // 384 cores on K620 * 5 = 1920

    // Initialize the type of thruster to be used
    thruster<double> thrust(0);

    //std::ofstream efficiencyGraph; // for viewing how many runge-kuttas ran per second for each combination of threads per block and total threads 
    //efficiencyGraph.open("efficiencyGraph.csv");
    std::cout << std::endl << "running optimize() with " << blockThreads << " threads per block and " << numThreads << " total threads";
    if (thrust.type) {
        std::cout << " using a NEXT-C thruster" << std::endl;
    }
    else {
        std::cout << " using no thruster" << std::endl;
    }

    optimize(numThreads, blockThreads, thrust); // optimize() currently declared in runge_kuttaCUDA.cuh
    //efficiencyGraph << blockThreads << "," << numThreads << "," << calcPerS  << "\n";
    //efficiencyGraph.close();
    
    delete launchCon;

    return 0;
}
