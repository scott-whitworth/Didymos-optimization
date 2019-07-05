"# MBS-optimization" 

**Scott Notes: 7/5/2019**

All files should have an over-arching header that describes what the file is doing. At least optimization.cpp should have project information identifying who is working on it and when.

Formatting consistancy: tabs should be 4 spaces. and indents should be consistant. If you don't know what I am talking about, I will show you. You can also use extensions to keep formatting consistant. 

optimization.cpp:
* Line 11/17: Why are we using timestamp()? Why are we using John Burkardt's? Is this notable for some reason?
* Line 23: Yes, I understand that type var = new type[n] is allocating memory, but what does start represent? What does step represent?
A better comment would be:  
`//Pre-allocating memory for starting parameters and steps
//start - starting parameters for optimization, order and contents defined in constants.h
//step - starting step sizes for starting parameters`  
* Line 28 - 74: Still think we should have a config file load in for the intial guesses, also probably other mission based parameters.
* Line 28 - 48: You can break these apart with some comments instead of relying on the offset keywords. Kind of readable, but also not at all obvious.
* Line 53 - 74: A big array with increasing index's and randome numbers. This is not great. At least keep the same formatting as the above parameters.
* line 51: 'initial change in variable size' - What does that mean? Are you sure that is what it is? What does changing this do?
* Probably want to throw some comments detailing the for loop at 76.
* Line 76: Why 10? That seems arbitrary
* Make sure you are deleting the correct start and step.
* Line 77: Are we going to want to vary step at all? It seems like if we are messing with start, we may also want to mess with step.
* Lines 85 - 89: This needs to be better and more descriptive. You are referencing things that are contained inside of start and step without ever defining start and step. Also fuel mass? Are you sure?
* Line 91: As we did not write the nelder mead.cpp, we should probably define these here as well, or at least our interpretation of them. i is not defined in nelder_mead.cpp
* line 102: This is dictating what this line is doing, not describing why we are doing it. a better comment would be `//Allocting xmin: space for nelder_mead algorithm to fill in final optimized parameters`
* Line 105: I think I know why this is here, but this is kind of what the whole program is doing. Maybe we would want `"Top of optimizing(): starting to calculate minimization of orbital motion"`
* Line 107-114: This is good information, but in the wrong spot and without context. Don't delete this but find a better location for this.
* Lines 118-123: Why are we picking these values?
* Lines 118-123: These would be better served to be assigned and defined in the inilization of these variables at line 95 - 98
* Lines 126 - 130 and Lines 142 - 146: This should be functionalized into a print parameters function. While we are at it, maybe we want to print out more than just a string of numbers? It would be far more useful to have well formatted output telling me what the 5th value means. This is also why these parameters really should be in a class (but that might not be clean with how NM requires an array). I think we could still overload the operator[] and get 'array' like operations without sacrificing the clarity of a class.
* Line 142: What is ifault? Is it an error code?
* Lines 155 - 164 should be functionalized ala trajectoryPrint()

acceleration.cpp/.h:
* This function could benefit from a mathematical overview of what is going on. What are the equations being used? That would go a long way to giving a clearer picture and solve a number of issues pointed out below.
* Comment in acceleration.h: massFuel is not the same as massExpelled? Also this calculates the magnitude of the accelleration, right? You define thrusting as being set in a different function, but what is it? 
* Line 6: For thrusting: Why is this a T? Could this be a flag? It is a little weird to compare a double to 0, that can cause issues if you are off by 0.00000000000000000000000000000001
* Line 3: thruster.h is already included in acceleration.h
* Line 9: I don't think this is true anymore. I might change Line 8 to be: `//If all of the fuel has been expelled, then no more thrust can be applied`
* Lines 14 and 15: This is not really relevant to the function. You have already defined(?) what thrusting is. You can leave it with `//No thrust will be calculated, coasting`
* Line 28: a better design would be to push this calcuation inside of the thruster. This looks true for all these lines, you need a buch of information from thrusterType. Might as well make this a method of thruster and pass in the relevant information (radius for the most part)
* Line 38: An example of a not-great comment. A better comment might be: `//update thrusterType's m_Dot based on power input`. This call does not return anything.
* Line 41: This comment is only part of the following line

calcFourier.cpp/.h
* Line 8 (.h): yes, by default this is an input array, what does it represent?
* Line 13 (.h): Probably should make the series a const ref. It kind of alreay is as an array, but still should probably make it a const at a minimum.
* Line 10/12 (.h) : If you are going to copy paste, make sure you copy and paste the right thing. This has nothing to do with gama or an in-plane angle. This just calculates the value of a series at a specific point in time.
Line 47 (.h): Seems like a good place for a bool. Be careful, this might change things, do this carefully.
* Lines 5-12: No comments. Should have at least a mathematical reference to what you are doing.
* Line 9: Slight optimization, it might be better to pre-calculate curTime/timeFinal. Depending on how often this is being called this could save a number of ops per series calculation. At the moment this needs to be computed twice * half the parameters per series calculation

constants.h
* Line 13: Torbital is not a great variable name for the orbtal period of the asteroid
* Line 14: Whoa... how is kiloConversion being used? Is there a case where this is not 0.001?
* Lines 17 - 34: No comments and lots of magic numbers. I get what some of these are, but they should be clearly documented, and these are not. Where did this data come from? It is in reference to a specific time right? What is that time?
* Lines 36 - 38: Official DART mission data that is two different values? That does not seem right. Additionally, is V_IMPACT a good name? That is the ideal V_IMPACT, not necessarily 'our' impact.
* Line 52: DRY_MASS naming consistancy.

elements.cpp/.h
* Lines 29 and 35 (.h): I don't think this is correct. Oh! This is more of an issue that these are not constructors, they are operator overloads.

motion_equations.cpp/.h

nelder_mead.cpp/.h

optimization.h

orbitalMotion.cpp/.h

runge_kutta.cpp/.h

thruster.cpp/.h






**Scott Notes: 6/19/2019**
* What is the difference between orbitalOptimization and development? These should probably be documented. They look like they are the same functionally, but significantly different.  
* Who is John Burkardt? (licensing?)
* optimizing is not the best name, and also why? This function should be better functionalized, this is not helped by the massive amounts of magic number setting at the top.
* every function needs a header, that header should be before the definition (declaration) not with the implementation
* What is nelder_mead? Why are we using it? Redundant to say 'see nelder.mead.cpp for i/o information. That should already be self-evident in your documentation. I lothe the use of 40 year old code.
* If we are using third party code, this should be clearly separated from our code. At the moment there is not really a demarcation.
* Why are we using a giant array? optimization lines 69-83? This is for sure going to lead to issues. There are much better ways to do this
* You have comments on these values, but they are not really defined in any meaningful way
* Our parameter space got really small, are we sure we are never going to go back and change? In general, we should optimized for a more difficult solution first, then pair back. Working from smaller numbers to bigger just punts the issue
* Bunch of magic numbers (69 - 109). These should not be in compiled code. These should be either rationalized, or loaded from a file.
* Formatting means something, we should used a standard format. Unless there is a specific reason, you should use 4 space 'tabs' for every indent.
* what is trajectory? Where is it defined?
* icount,numres and ifault are never initialized
* `"\n"` can be in the same `cout << "\nF(X*)" << ynewlo ` line  
 * Where is trajectoryPrint defined?  

 nelder_mead.cpp:  
 * This does not look like the original source code, if we modified it, we need to make sure we know where. If someone else modified it, we need to make sure we are referencing that. This code suggests it is modified FORTRAN from the 70's. C++ style new dynamic memory was not introduced until ~1983.
 * Additionally, terrible commenting.
 * I am pretty sure at the end of the day we are not going to use this, but what is important about this is not *if* it works, but *how* it works. This is not a code issue, but a project path issue. 
 * This is a great example of why functionalization is so important. If you ever want to make a change, it is nearly impossible due to the complexity of the nested loops/ifs
 * Also a great example of terrible variable names. The line ` l = l + 1;` is not only non-descriptive, it is really easy to look at as `1 = 1 + 1;`
* Why are we using their timestamp? Didn't we already write a time stamping function? Ugg... any why are we using c-style strings.

 orbitalMotion.h/.cpp: 

* trajectoryPrint seems like it is doing the same thing as trajectory. This will cause issues. Better functionalization is needed.
* If we are calling trajectory() a bunch from within the nelder_mead function, we should make sure this is free of any allocations and printouts, these take time.
* For readability: if you have a variable/function that is not defined in a function, it should be obvious where that variable/function is coming from. orbitalMotion.cpp:32. This can be in a comment or a declaration at the top of the file
* Line 39 is a great example of where functionalization would help
* 68-74 woof. This is an intense way to do this. Not at all obvious or documented.
* It seems like absTol and Fmin are constants. They should be declared as such.
* *Possible major error* I think you are checking the cost function twice. Once in the trajectory() function and once in the nelder_mead function. Both are being checked against different values. Are you sure you want to be doing that?
* Optimization opportunity: earthInitial will calculate similar values many many times. There is probably a better way for us to keep track of values we have already calculated.
* What is the difference between rk4Simple and rk4sys? If there are two funcions, there is another vector for issues.

Constraints:
* These constraints should be more obvious. If we are going to be using a #define for an offset, this needs to be in the form of TRIPTIME_OFFSET. It is not obvious at orbitalMotion.cpp:29 that tripTime_offset is not 'just' another variable, but instead a const offset.







