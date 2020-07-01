<h1>DART Mission Optimization Project</h1>
Last updated: July 1st, 2020

<h2>Project Background & Current Objective</h2>

  NASA's Double Asteroid Redirection Test (DART) mission involves having the spacecraft perform a kinetic impact to change the orbital trajectory of Dimorphos around the parent body, Didymos.  The spacecraft to be used is fitted with a NEXT ion thruster though is not required to hit the target but rather be used as a technical demonstration.  For this project, our goal is to find an optimal trajectory using the thruster that would result in a greater difference in velocity with the target to make for a more effective change in Dimorphos' orbit.  Along with optimizing the effectiveness of the DART mission, this project also aims at optimizing the process of finding a solution.

  Currently, the tools for finding the best trajectory is utilizing a genetic algorithm that uses Nvidia's CUDA platform to optimize parameters that leads to both hitting the asteriod and maximizing the velocity difference.  At this stage of development the focus is on improving and verifying the effectiveness of the genetic algorithm in finding valid parameters that leads to an impact on the asteriod.

  Parameters being optimizing are the following
  | Variable Name               | Units    	  | Description                                                                                                                                                |   	|
  |----------------------------	|------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
  | Trip Time                  	| Seconds (s) | How long of a duration the trip for the spacecraft takes from once it leaves Earth's sphere of influence to impact with Dimorphos, the impact date is predetermined by NASA but launch date is flexible and so this variable also impacts the initial position and velocity of the craft when it leaves Earth's sphere of influence|   	|
  | Alpha launch angle          | Radians  	  | The in-plane angle of the spacecraft as it leaves Earth's sphere of influence                                           |   	|
  | Beta launch angle           | Radians  	  | The in-plane angle from tangent of Earth's sphere of influence that the spacecraft's intial velocity takes |   	|
  | Zeta launch angle           | Radians  	  | The out-of-plane angle of the spacecraft's initial velocity when leaving Earth's sphere of influence |   	|
  | Gamma coefficients          | None  	    | Used in a fourier series to calculate the gamma (in-plane angle in radians) angle at a given time for thruster behavior |   	|
  | Tau coefficients            | None  	    | Used in a fourier series to calculate the tau (out-of-plane angle in radians) angle at a given time for thruster behavior                                           |   	|
  | Coast coefficients          | None  	    | Determines if the thruster is activated or not (coasting) at a given time for thruster behavior |   	|

<h2>Files & Folders Navigation Guide</h2>
There are many folders and files in this project, so here's a short rundown of folder contents to help navigate everything;

  - Cuda: Where the most recent optimization code that attempts to find a best trajectory can be found, uses the CUDA platform to use  GPU and genetic algorithm 
    * Config_Constants: Where cudaConstants structure is defined and default genetic.config file is, cudaConstants handles storing const values that we may want to be able to change for different runs of the program.  Dependent on constants.h file
    * Earth_calculations: Code for calculating the earth conditions and defines the global pointer variable launchCon (earthInfo.h). Dependent on Motion_Eqns/elements.h, Thrust_Files/thruster.h, and Config_Constants/config.h.
    * Genetic_Algorithm: Defines individuals used in the genetic algorithm and crossover methods to generate new generations in a pool.  Contains main file (titled optimization.cu) that has main() and main optimize function that is called.
    * Motion_Eqns: Defines elements structure that is used to describe the position and velocity of an object in space (such as Earth). Dependent on Thrust_Files and Config_Constants/config.h.
    * Runge_Kutta: Holds the runge_kutta functions with versions for both CPU and GPU usage. Also defines rkParameters structure that holds the information that the genetic algorithm attempts to optimize.  Dependent on Thrust_Files, Config_Constants, etc.
    * Thrust_Files: Contains the code that describes the thruster and its behavior using a Fourier series to determine angles of direction. Dependent on Config_Constants/config.h
    * constants.h: Stores constant properties, such as AU unit value and optimized variable offsets for the array that stores the values, these are constants that should not be easily changed.
    * optimizedVector.bin: Binary file containing 14 parameters that can be used as initial guesses in the genetic algorithm, values derived from previous code in orbitalOptimization folder which is also based on old impact date data.
  
  - PostProcessing: Contains MATLAB files that take in output files from the Cuda program to display results.
  - archives: Contains MATLAB code files that are considered no longer relevant to the project but shouldn't be completely discarded.
  - cuda_framework: Contains code that was used in setting up and testing usage of the CUDA platform, not actively used in finding a solution
  - earthPosition: Isolated version of earthInfo used in Cuda to calcualte the path of the Earth, not updated to usage of cudaConstants config file reading
  - orbitalOptimization: 2019's non-GPU dependent optimization method that utilizes Nelder-Mead instead of a genetic algorithm, also not actively used in current appraoch to the task.
  - DerivingImpactValues.xlsx : A current (likely temporary) method of converting JPL data values from cartesian to cylindrical units by using spreadhsheet fields and computations.  The converted values are then copied and put into genetic.config to be read and used.  This does not directly impact the code.

<h2> Running CUDA Code: </h2>

On WU System & Personal Computers:
1. To Compile:
   1. Open VSCode and open the Didymos Optimization project folder
   2. Open command prompt in the VScode terminal.
   3. Navigate to the Genetic_Algorithm folder (input and enter "cd Cuda" then "cd Genetic_Algorithm").
   4. Enter the following
      `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.22.27905\bin\Hostx64\x64\cl.exe" -o <OutputFileNameHere>.exe optimization.cu`   
    - Make sure that <OutputFileNameHere> has a specific name (it will create a .exe file with whatever you named it), a file path could also be made so that the files created from the compiler are not in the same folder
    - *If your computer can’t find the file, then there is a problem with the path. Copy and paste the file path into a file explorer, and delete everything after `\MSVC\` then click through until you get to `\cl.exe`. Copy this path and use it to replace the old path.*
    - This code should add an .exe file with the output name that was entered

  For Tesla Machine:			  
   - Do the same as above, but for step 3 use this command:
     `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl.exe" optimization.cu -o <nameOfOutput>`
		
2. To Run:
    1. Open command prompt in the VScode terminal.
    2. Navigate to the Genetic_Algorithm folder (input and enter "cd Cuda" then "cd Genetic_Algorithm").
    3. Type the .exe file name from above (don't forget to add .exe) and enter
    4. The program will begin and should show the following in order on the terminal;
       1. Outputs the GPU device name and intial values read from genetic.config that is in Config_Constants folder.
       2. Calculate the Earth data with a visible loading bar.  The range of this data is 0.5-2.5 years before impact date
       3. Outputs the number of threads and blocks that will be used in the optimize function and starts the algorithm.
       4. On the terminal, displays a "." for every generation calculated and sorted.  Every disp_freq generation it displays the current generation number (how many have been calculated up to this point minus 1) and best individual in the pool.  Also displays change in anneal size every change_check generations.
       5. Along with terminal display, there are serveral file outputs made during the program's run.
       6. Once the best_count number of individuals in the pool have passed the tolerance (less than pos_threshold set in config file), the algorithm has "succeeded" and will output best_count number of files that describe that individual's parameters that can be used in PostProcessing to observe.
       7. The program is finished and so closes

3. Changing properties:
      In Config_Constants, genetic.config holds a set of variables that can be changed before running the .exe file.  Refer to config_readme.md for specifics on how each variable impacts the program's behavior and format of the file.

4. Output Resulting Files
  - EarthCheckValues.csv : Recorded prior to starting the optimizing genetic algorithm, contains the position and velocity for every day of Earth (lower resolution than what is calculated which is every hour).  Used to verify that it matches with the JPL database.
  - BestInGenerations.csv/.bin : A record of the best individual's parameters (not including thruster values) and posDiff and velDiff as well as the following generations anneal value.  Records every write_freq generation.  Also records the best individual at the end of the algorithm.
  - WorstInGenerations.csv/.bin : Similar to BestInGenerations except it selects the worst individual instead.
  - final-optimization#.bin and orbitalMotion-accel#.bi : For best_count number of individuals, files are made of this format when the algorithm is finished to save the best individuals (replace # with and you get that numbered best individual) to then be used by PostProcessing to show the best trajectory found.

<h2> Running Matlab Code </h2>

0. In order to run the Matlab code, you must be on a Whitworth machine that has Matlab installed.
1. Copy the binary files outputted by the Cuda application generated when it finished.
2. Open Matlab client with directory to PostProcessinng contents.
3. In the command window, enter "filePlot(#)" with # as the number referring to the individual being plotted (1 would be the best individual, 2 would be 2nd best, etc.).
4. Graphs will be generated that show the path in a three dimensional graph, coast behavior, etc. that could be exported into png format.

<h2>NASA JPL Data for Impact Date Position & Velocity</h2>
Here is how the impact data was obtained to be used in the config value.

1.  Navigate to https://ssd.jpl.nasa.gov/horizons.cgi, this database is used to retrieve final position and velocity vector components for both Earth and Didymos relative to the Sun.

2.  The ephemeris type should be a vector table. In table settings, select Type 2. Set the coordinate origin to the Sun's center of mass. The current impact date is 09-30-2022 19:54:55 UTC, so set an adequate time window with a resolution of 1 minute. Set the units to a.u. for position and a.u./day for velocity. Then, select Generate Ephemeris.

3.  To obtain final Earth elements, change the target body to Earth-Moon barycenter and generate a new ephemeris.

4.  Using impactParams.m, DerivingImpactValues.xlsx, or some other equivalent method, convert the values to cylindrical coordinates with velocity values changed from AU/day to AU/s.