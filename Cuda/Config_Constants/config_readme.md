<h1> Config File Specifications/Information </h1>
Last Updated: July 13th, 2020

<h2>File format for config file</h2>

- The config file allows empty rows and comments ("//" at start of comment line) for formatting the presentation of the contents, currently does NOT allow in-line comments or spaces in variable assignments
- The parsing process currently does not attempt any verification of assigning values to variables (lack of assignment nor duplications).
- When reading the file, the assumption is made that the config file contains valid values for all variables and will result in exception thrown if invalid value (such as string instead of double) is used

<h2>The cudaConstants struct</h2>

- In the code, the structure that uses the config file is called <b>cudaConstants</b> and is only accessed when being constructed (therefore changing the config file during a run would have no impact).
- An overloaded << operator for cudaConstants that outputs the object's contents with labelling/formatting for better readibility to outputt onto terminal screen in the main optimization program.
- Function compareConstants() takes in two const cudaConstants that returns true if all variables are equivalent, used in the genetic algorithm as means of verifying that the values are not changing during runtime.
- For changing what config file is used, the file address can be changed where cudaConstants is declared within the main function in optimization.cu.
- Default address is "genetic.config" in same folder as the .exe file, optimization.cu has address set as "../Config_Constants/genetic.config".
- If config file address is invalid, will output to terminal that is the case.

<h2>Variables in Config/cudaConstants</h2>

Table 1. Setup & General Values
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| time_seed                  	| int/string 	| None  	| Sets the seed used in random generation, either specify a seed to use or place "NONE" for the seed to be time(0)                                           	                    |   	|
| rng                     	    | mt19937_64    | None  	| Random number generator that uses time_seed as initialization value                                                          	                                                    |   	|
| random_start               	| boolean    	| None  	| If "true", sets initial generation's individuals to hold parameters with random values, if "false" it initializes the individuals from a provided file                  	        |   	|
| initial_start_file_address 	| string     	| None  	| If random_start is false, the program uses this address to get parameter values for the initial individuals with the assumption that the file hold 14 sets 	                    |   	|
| pos_threshold              	| double     	| AU      	| Sets the maximum positional difference of the spacecraft to the target at end of its trajectory path                                                    	                        |   	|
| write_freq                 	| int        	| None  	| Sets number of generations to process before writing information onto files, 1 is to write every generation                                               	                    |   	|
| record_mode                 	| boolean       | None  	| If "true", sets program to output various files that describe the performance, meant to be used in helping verify/debug behavior.  Currently some file output does not properly support non-default parameter sizes |   	|
| disp_freq                  	| int        	| None  	| Sets number of gnerations to process before outputting to console terminal, 1 is to display output every generation                                       	                    |   	|
| change_check               	| int        	| None  	| For how many generations until it checks to see if the best individual has changed, if no change the anneal value is reduced by multiplying with anneal_factor                    |   	|
| rk_tol                 	    | double     	| None  	| The relative/absolute (not sure which one it is) tolerance for the runge kutta algorithm	                                                                                        |   	|
| f_min                 	    | double     	| None  	| The expected precision for the optimization cost convergance. This number is meant to avoid unnecesary iteration whitin neder _ mead	                                            |   	|
| max_numsteps                 	| double     	| None  	| Used for time stepping in runge_kuttaCuda.cu	                                                                                                                                    |   	|
| startTime                 	| int        	| seconds   | The time offset backwards from impact date values of Earth-Moon Barycenter to start reverse runge-kutta calculating positiona and velocity	                                    |   	|
| durationTime                 	| int        	| seconds   | The time duration backwards from the offset of Earth-Moon Barycenter to start reverse runge-kutta calculating positiona and velocity, used in deriving endTime	                |   	|
| endTime                    	| int        	| seconds   | The time offset backwards from the impact date of Earth-Moon Barycenter to end reverse runge-kutta calculating positiona and velocity, derived from (startTime + durationTime)    |   	|
| timeRes                    	| int        	| seconds   | The "gap" between each calculation for Earth's backward runge-kutta, for example 3600 sets every calculation to be 1 hour apart                                                   |   	|

Table 2. Genetic Algorithm Values
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| num_individuals           	| int        	| None  	| Sets the size of the population pool and number of threads used as an individual is given a thread, recommended to not change 	                                                |   	|
| survivor_count               	| int        	| None  	| Number of individuals selected as "survivors" to produce new individuals in the next generation in the genetic algorithm, every pair produces 8 new individuals, value must be even|   	|
| thread_block_size           	| int        	| None  	| Number of threads per block on the GPU being used, recommended to not change 	                                                                                                    |   	|
| anneal_initial             	| double     	| None  	| The initial anneal value used, anneal impacts the maximum possible mutation value when generating a new individual (does not impact probability) 	                                |   	|
| anneal_factor             	| double     	| None  	| The multiplier applied to anneal value if no change in the best individual is occurring                                                                        	                |   	|
| best_count                 	| int        	| None  	| How many individuals must have obtained a solution before ending the algorithm, also outputs the top number of individuals up to best_count 	                                    |   	|
| mutation_rate              	| double     	| None  	| The probability of a mutation occurring when generating a new individual, gurantees at least one gene is changed                                           	                    |   	|
| double_mutation_rate       	| double     	| None  	| Probability that if a mutation is occurring that it affects two genes 	                                                                                                        |   	|
| triple_mutation_rate       	| double     	| None  	| Probability that if a mutation is occurring that it affects 3 genes 	                                                                                                            |   	|
| gamma_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for gamma values 	                                                                                                                            |   	|
| tau_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for tau values 	                                                                                                                            |   	|
| coast_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for coast values 	                                                                                                                            |   	|
| triptime_mutate_scale 	    | double     	| None  	| Affects the maximum mutation range for triptime values 	                                                                                                                        |   	|
| zeta_mutate_scale          	| double     	| None  	| Affects the maximum mutation range for zeta values 	                                                                                                                            |   	|
| alpha_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for alpha values 	                                                                                                                            |   	|


Table 3. Mission Values
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| thruster_type                	| int        	| None  	| Determine what thruster is used, 0 for none and 1 for NEXT ion thruster 	                                                                                                        |   	|
| dry_mass                     	| double        | kg      	| Set the mass of the spacecraft without fuel, also used in determining wet_mass 	                                                                                                |   	|
| fuel_mass                     | double        | kg      	| Sets the initial mass of fuel in the spacecraft, used in determining wet_mass 	                                                                                                |   	|
| wet_mass                     	| double        | kg      	| The total mass of the spacecraft with fuel                                                                                                                  	                    |   	|
| coast_threshold             	| double     	| None  	| In a range from 0 to 1, 1 sets the spacecraft to coast at all times while 0 sets the spacecraft to always have thruster on 	                                                    |   	|
| c3energy                     	| double     	| m<sup>2</sup>/s<sup>2</sup>  	| The specific energy of the spacecraft when leaving the sphere of influence of the earth-moon center of mass, determines the magnitude of the escape velocity that is stored in v_escape 	|   	|
| v_escape                     	| double     	| AU/s  	| The magnitude of the initial velocity of the spacecraft when leaving the sphere of influence of the earth-moon center of mass, not in config file but rather derived from c3energy 	                        |   	|
| v_impact                 	    | double     	| AU/s  	| NASA's official mission impact velocity difference the spacecraft will collide with Dimorphos, does not impact the performance of the code	                                    |   	|


Table 3a. Impact Position & Velocity Values
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| r_fin_ast           	        | double     	| AU      	| The radius position of the target at impact date, relative to the Sun bodycenter 	                                                                                                        |   	|
| theta_fin_ast        	        | double     	| Radians  	| The theta angle position of the target at impact date, relative to the Sun bodycenter	                                                                                                    |   	|
| z_fin_ast           	        | double     	| AU      	| The z (off-plane offset) position of the target at impact date, relative to the Sun bodycenter	                                                                                        |   	|
| vr_fin_ast           	        | double     	| AU/s  	| The velocity of the radius component of the target at impact date, relative to Sun bodycenter	                                                                                        |   	|
| vtheta_fin_ast      	        | double     	| AU/s  	| The tangental velocity of the target at impact date	                                                                                |   	|
| vz_fin_ast           	        | double     	| AU/s  	| The velocity of the z component of the target at impact date, relative to the Sun bodycenter	                                                                                            |   	|
| r_fin_earth           	    | double     	| AU     	| The radius position of the earth-moon center of mass at impact date, relative to the Sun bodycenter	                                                                                                            |   	|
| theta_fin_earth          	    | double     	| Radians  	| The theta angle position of the earth-moon center of mass at impact date, relative to the Sun bodycenter	                                                                                                    |   	|
| z_fin_earth           	    | double     	| AU      	| The z (off-plane offset) position of the earth-moon center of mass at impact date, relative to the Sun and used to plot it's path backwards in time for launch positions of the spacecraft 	        |   	|
| vr_fin_earth           	    | double     	| AU/s  	| The velocity of the radius component of the earth at impact date 	    |   	|
| vtheta_fin_earth         	    | double     	| AU/s  	| The tangental velocity of the earth-moon center of mass at impact date 	|   	|
| vz_fin_earth           	    | double     	| AU/s  	| The velocity of the z component of the earth-moon center of mass at impact date, relative to Sun bodycenter and used to plot it's path backwards in time for launch positions of the spacecraft 	            |   	|
