<h1>DART Mission Optimization Project</h1>
Last updated: June 23rd, 2020

<h2>Project Background & Current Objective</h2>
  NASA's Double Asteroid Redirection Test (DART) mission involves having the spacecraft perform a kinetic impact to change the orbital trajectory of Dimorphos around the parent body, Didymos.  The spacecraft to be used is fitted with a NEXT ion thruster though is not required to hit the target but rather be used as a technical demonstration.  For this project, our goal is to find an optimal trajectory using the thruster that would result in a greater difference in velocity with the target to make for a more effective change in Dimorphos' orbit.  Along with optimizing the effectiveness of the DART mission, this project also aims at optimizing the process of finding a solution.

  Currently, the tools used for finding the best trajectory is utilizing a genetic algorithm that uses Nvidia's CUDA platform to optimize parameters that leads to both hitting the asteriod and maximizing the velocity difference.  At this stage of development the focus is on improving and verifying the effectiveness of the genetic algorithm in finding valid parameters that leads to an impact on the asteriod.

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

<h2>Current Approach</h2>
--WIP-
  The program utilizes a config file to determine multiple properties that affect the performance of the algorithm.  This includes but not limmited to mutation rates and how often to record the progress of its best individual.  It also includes options such as the asteriod's location and mass of the spacecraft.  For more information on the config file, refer to <b>config_readme.md</b> found in CUDA/Config_Constants.
  The genetic algorithm as of the time of this writing does not optimize the trajectory for higher speed velocity, currently only concerned with reducing the positional difference to what is specified in the config file under pos_threshold.

<h2> Running CUDA Code: </h2>

On WU System:
1. To Compile:
   - Enter the following into command prompt within VSCode 

      `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.22.27905\bin\Hostx64\x64\cl.exe" -o <OutputFileNameHere>.exe optimization.cu`

   
   - Make sure that <OutputFileNameHere> has a specific name (it will create a .exe file with whatever you named it)

    
     *If your computer can’t find the file, then there is a problem with the path. Copy and paste the file path into a file explorer, and delete everything after `\MSVC\` then click through until you get to `\cl.exe`. Copy this path and use it to replace the old path.*

    - This code should add an .exe file with the output name that was entered
		
  1. To run:
			
      Type the .exe file name from above (don't forget to add .exe) and enter
      
  2. Changing properties:
    
      In Config_Constants, genetic.config holds a set of variables that can be changed before running the .exe file.  Refer to config_readme.md for specifics on how each variable impacts the program's behavior.
      
		
  For Tesla Machine:
			  
   - Do the same as above, but compile using this:

     `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl.exe" optimization.cu -o <nameOfOutput>`
       
 
 <n>
<h2> Running Matlab Code </h2>

  