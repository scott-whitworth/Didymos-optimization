 <h1> Asteroid Impact Mission 

<h2> Running CUDA Code: </h2>

 On WU System:
1. To Compile:
   - Enter the following into command prompt within VSCode 

      `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.22.27905\bin\Hostx64\x64\cl.exe" -o <OutputFileNameHere> optimization.cu`

   
   - Make sure that <OutputFileNameHere> has a specific name (it will create a .exe file with whatever you named it)

    
     *If your computer can’t find the file, then there is a problem with the path. Copy and paste the file path into a file explorer, and delete everything after `\MSVC\` then click through until you get to `\cl.exe`. Copy this path and use it to replace the old path.*

    - This code should add an .exe file with the output name that was entered
		
  2. To run:
			
      Type the .exe file name from above (don't forget to add .exe) and enter
      
  3. Changing properties:
    
      In Config_Constants, genetic.config holds a set of variables that can be changed before running the .exe file.  Refer to config_readme.md for specifics on how each variable impacts the program's behavior.
      
		
  For Tesla Machine:
			  
   - Do the same as above, but compile using this:

     `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl.exe" optimization.cu -o <nameOfOutput>`
       
 
 <n>
<h2> Running Matlab Code </h2>

  





