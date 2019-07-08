"# MBS-optimization" 

**Scott Notes: 7/5/2019**

optimization.cpp:
* Line 28 - 74: Still think we should have a config file load in for the intial guesses, also probably other mission based parameters.
* Lines 118-123: These would be better served to be assigned and defined in the inilization of these variables at line 95 - 98
* Lines 126 - 130 and Lines 142 - 146: This should be functionalized into a print parameters function. While we are at it, maybe we want to print out more than just a string of numbers? It would be far more useful to have well formatted output telling me what the 5th value means. This is also why these parameters really should be in a class (but that might not be clean with how NM requires an array). I think we could still overload the operator[] and get 'array' like operations without sacrificing the clarity of a class.










