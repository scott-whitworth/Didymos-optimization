<h1> Config File Specifications/Information </h1>
(Current note from June 22nd) Needs to be updated with new variables used in program

<h2>The config file allows empty rows and comments ("//" at start of comment line) for formatting the presentation of the contents</h2>
<h2>For changing what config file is used, the file address can be changed where gConstant is declared within the main function</h2>
<h2>If config file address is invalid, will output to terminal that is the case.  Also assumption is made that the config file contains valid values for all variables</h2>
<h2>As of writing (June 17th, 2020), the program is set to require the file genetic.config to be in the same folder as the executable and additional note that the optimizeVector.bin file must be one folder outside of the executable's (address is "../optimizedVector.bin") </h2>

<h2>Future Considerations</h2>
Consider expanding on either current structure/config to contain elements beyond the genetic algorithm (mission parameters for example) or have multiple config files for the various components.  Also expansion to have the program allow easy changes to other config files to allow multiple versions without needing to recompile with address changes (user prompt?  A config file for choosing config files?).  Consider having the address for optimizedVector.bin either be within same folder as executable, or have these dependencies (or where file outputs are stored) also be malleable variables within the config files.


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