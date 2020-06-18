<h1> Config File Specifications/Information </h1>

<h2>The config file allows empty rows and comments ("//" at start of comment line) for formatting the presentation of the contents</h2>
<h2>For changing what config file is used, the file address can be changed where gConstant is declared within the main function</h2>
<h2>If config file address is invalid, will output to terminal that is the case.  Also assumption is made that the config file contains valid values for all variables</h2>
<h2>As of writing (June 18th, 2020), the program is set to require the file genetic.config to be in the folder Config_Constants that is a neighboring folder to where the the executable is located (address is "../Config_Constants/genetic.config") and additional note that the optimizeVector.bin file must be one folder outside of the executable's (address is "../optimizedVector.bin") </h2>

<h2>Future Considerations</h2>
Consider expanding on the config holding more variables that impact the algorithm or its conditions.  Also expansion to have the program allow easy changes to other config files to allow multiple versions without needing to recompile with address changes (user prompt?  A config file for choosing config files?).  Consider having the address for optimizedVector.bin either be within same folder as executable, or have these dependencies (or where file outputs are stored) also be malleable variables within the config files.


<h2>Variables & Description</h2>

- time_seed: Sets the randomization seed used in optimize(), if desiring truly random replace numeric value with NONE (will output onto terminal if seed is set to time(0))

- write_freq: Sets for how gap between number of generations before the best individual's properties are appended to BestInGenerations.csv and BestInGenerations.bin (the last generation is always appended at the end of the file)
- disp_freq: Sets for how gap between number of gnerations before outputting information to the terminal, does not impact the algorithm's process of finding a solution or recording the progress 

- best_count: How many individuals are required to reach the criteria of a solution before ending the program and recording the solutions into binary files
- pos_threshold: Sets the maximum positional difference of a solution
- speed_threshold: Sets the minimum speed difference of a solution
- change_check: How often the algorithm checks to see if the best individual has changed, if no change then the anneal value is changed
- anneal_factor: Multiplier for changing the anneal value if no change in best individual
- anneal_initial: The initial anneal value for the algorithm


- mutation_rate: The probability percentage of mutations occurring
- double_mutate_rate: The probability percentage of a mutation affecting 2 genes in total
- triple_mutate_rate: The probability percentage of a mutation affecting 3 genes in total
- gamma_mutate_scale: Affect the max range for mutating the gamma values
- tau_mutate_scale:: Affect the max range for mutating the tau values
- coast_mutate_scale: Affect the max range for mutating the coast values
- triptime_mutate_scale: Affect the max range for mutating the triptime value
- zeta_mutate_scale: Affect the max range for mutating the zeta value
- beta_mutate_scale: Affect the max range for mutating the beta value
- alpha_mutate_scale: Affect the max range for mutating the alpha values

- thruster_type: Chooses what thruster is used
- c3energy: The specific energy of the spacecraft when exiting the earth's sphere of influence, determines the maganitude of the initial velocity

- earth conditions (such as r_fin_earth, vtheta_fin_earth): The final position and velocity of earth on impact date, used in determining earth's positions across the possible launch window time-frame
- asteriod conditions (such as r_fin_ast, vz_fin_ast): The final position and velocity of the asteriod on impact date, used in determing the target location and velocity to compare with results of the spacecraft's trajectory