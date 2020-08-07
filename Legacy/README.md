# Legacy Code for Reference

*Should not be used, only for reference.*

## Summer 2019 Work
Work developed by Lauren, Ben and Mateo
* `orbitalOptimization` - CPU implementation of exiting Matlab code. 
  * Utilizes Runge Kutta algorithm based on Matlabs ODE45 library
  * Nelder Mead method for minimizing final position difference
  * Heavily developed to test the implementation in C++, achieved significant performance improvement
* `earthPosition` - Development of Earth Initial position calculation for use with the GPU code
  * One optimized variable is trip time
    * Starting position is computed backwards in time from the impact date
  * Each step / individual has a different trip time, thus many different initial earth positions are needed
  * Integrated into the CUDA code as earthInfo.h/.cpp
* `cuda_framework` - first development of GPU based RK calculations.
  * very basic sandbox for setting up kernel and GPU based backend
  * eventually formed the first optimization.cu file that included genetic algorithm 
* `archives` - early work on code development



