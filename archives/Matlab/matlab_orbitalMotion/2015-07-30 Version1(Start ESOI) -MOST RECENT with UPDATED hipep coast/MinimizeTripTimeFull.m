%MINIMIZETRIPTIME Searches for the lowest converging trip time guess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is designed to take an order for a fourier
%fit, an estimate of how long the trip will take and then a convergance
%radius for its solution. Given that information it will optimize
%trajectories for the ship to land on the asteroid for successively smaller
%time guesses until it finds the minimum possible time in which the ship
%can successfully land. That value is then returned as the most efficient
%time. Much of the relevant problem information specifically about the
%ship's pecifications and asteroid orbital parameters are read in from
%other files. Almost all functionality of the program branches off of this
%function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MultipleRuns
function [finalTime,escapeTime,timeToExecute,phiCoeff,plotFileName,rError,thetaError,uError,hError,fval,lastConverge,massPowerPlant,massThruster,massStruct,mdot,lowestEscape,lowestConverge,lowestReturnConverge,lowestForwardConverge,refinedTimeStep,FuelMass,lowestdepartBound,lowestapproachBound,lowestcoastFraction] = MinimizeTripTimeFull(FS,order,direction,massFuel,numberOfEngines,alphaP,alphaT,inputPower,uExhaust,efficiency,massPayload,massStruct,massSample,lastConverge,selectForwardAltitude,selectReturnAltitude,entryVelocity,entryAngle,astDes,timeWait,isSolar,powerLevel,thrusterName,numPoints,thrusterStructs,maxTimeLimit,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,flagStore,optimize,interpolateStruct,throttleFlag)

%starts timer for function
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non-dimensionalize and initialize constants

%Non-dimensionalize constants
earthRadius = Constants.EARTHR / Constants.MCONVERSION;
forwardAltitude = selectForwardAltitude / Constants.MCONVERSION;
returnAltitude = selectReturnAltitude / Constants.MCONVERSION;
entryVelocity = entryVelocity * (Constants.SCONVERSION/Constants.MCONVERSION);
entryAngle = entryAngle * Constants.DEGTORADCONVERSION;
%Initialize variables

%Non-dimensional starting altitude of launch measured from earth's center
forwardRadius = earthRadius + forwardAltitude;

%Non-dimensional ending altitude of spacecraft measured from earth's center
returnRadius = earthRadius + returnAltitude;

%Initialize mass of the power plant and the thrusters
massPowerPlant = numberOfEngines * alphaP * inputPower;
massThruster = numberOfEngines * alphaT * inputPower;

%Initializes final values of r, theta, and u for Earth and Asteroid, and h
%of Asteroid based on text files stating the relative positions of earth
%and asteroid. It calculates the point of closest approach and sets the
%conditions of closest approach as the conditions at landing.
[earthCloseApp,asteroidCloseApp,AsteroidSemimajor,AsteroidEccentricity] = GetFinalState(astDes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The forward trip also carries the fuel for the return trip as calculated
%in MultipleRuns.
massPayload = massPayload + massFuel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fminsearch is an unconstrained optamization package in MATLAB.  It gets
%passed a vector of independent variables, in this case 1 variable
%(timeGuess) and a function to evaluate each guess with. The
%maximum number of function evaluations and maximum iterations must not be
%set too high because fminsearch has a tendency to get trapped at random
%time guesses and will iterate infinitely. In almost every case, fminsearch
%will produce an F value in less than 20,000 iterations, or it will iterate
%infinitely. The tolerance on the independent variable (TolX) and the
%tolerance on the function evaluation value (TolFun) in this circumstance
%will be very close to identical.  If the ship successfully lands on the
%asteroid in the time allotted, then timeObjectiveFunction will return the
%timeGuess as the function value. If it fails, timeObjectiveFunction
%will return a very large value.  As long as the ship is landing on the
%asteroid, the variation in X and the objective function will be identical,
%Therefore, they are assigned to the same value.

%A drawback of using fminsearch like this on a 1 dimensional function is
%that it knows nothing about the search space.  We know that the next
%lowest answer will always be at a lower time guess than the previous but
%there is no way to input that information to the algorithm. For some time
%we used a binary search but from time to time there are irregularities in
%the results from timeObjectiveFunction and fminsearch is much more capable
%of handling those than our binary search was.

%As a precurser to the main time optimizing search algorhythm, we wrote an
%algorithm to generate a good starting time guess from which the main
%algorithm begins.  This is accomplished by beginning at the lowest
%feasible time that a trip could take and stepping at a large time interval
%(~0.05 yr) forward through time until the value returned by fminsearch for
%the value of the function minimum dips below the 'IsInteresting' threshold
%(1E-10). Once such a value is found, the algorithm steps backward at a
%smaller time interval (~0.01 yr) until it finds two instances of values
%for the function returned by fminsearch which are above the threshold. The
%time at which this occurs is stored as the value for which the main
%time optimization algorithm will begin.  From there, the algorithm steps
%forward in time an even smaller time step until it reaches a converged
%value (1E-20).  From there, the algorithm will take progressively smaller
%backwards steps to try and find a lower converged value.

%We also wrote a new binary search MinimizeTripTime as well, which double checks
%around mostly all non-converged points to assure that fminsearch has found 
%the best possible optimization. The main goal of this was to decrease
%computation time by eliminating the need for unnecessary time guesses.

lastfval = nan;
lastEscape = nan;

%Setting the max and min time guesses
if numPoints == 1
    largeForwardStep = rand(1) * .01 + 0.0375; %interval for large forward steps
    smallerBackStep = -0.01; %interval for smaller backward steps
    if direction == Direction.RETURN
        lowerTimeLimit = rand(1) * .1 + .6; %lowest expected limit of trip time
    else
        lowerTimeLimit = rand(1) * .1 + .6; %lowest expected limit of trip time
    end
    upperTimeLimit = 20; %highest expected limit of trip time
    innerLowerTimeLimit = lowerTimeLimit;
    if direction == Direction.RETURN
        lastConverge = nan;
        maxTimeGuess = 20;
        timeStep = .007 + rand(1)*.001;
    elseif direction == Direction.FORWARD
        maxTimeGuess = 20;
        timeStep = .007 + rand(1)*.001;
    end
else
    newStartFraction = 0.75;
    largeForwardStep = ((1-newStartFraction)* maxTimeLimit)/5; %interval for large forward steps-making it proportional to the previous trip time
    smallerBackStep = -largeForwardStep/2; %interval for smaller backward steps
    lowerTimeLimit = maxTimeLimit * newStartFraction; %Lowest expected trip time
    upperTimeLimit = 20;
    innerLowerTimeLimit = .01 + rand(1)*.001;
    if direction == Direction.RETURN
        lastConverge = nan;
        maxTimeGuess = 20;
        timeStep = -smallerBackStep/5 + rand(1)*(-smallerBackStep/10);
    elseif direction == Direction.FORWARD
        maxTimeGuess = 20;
        timeStep = -smallerBackStep/5 + rand(1)*(-smallerBackStep/10);
    end
end

%Initialize lowest values
lowestConverge = maxTimeGuess;
lowestcConverge = (pi/2)*(-1 + 2*rand(1,2*order+3));
lowestfval = maxTimeGuess;
lowestfinalMass = maxTimeGuess;
lowestescapeVelocity = maxTimeGuess;
lowestescapeEarthAngle = maxTimeGuess;
lowestearthConditions = zeros(4,1);
lowestasteroidConditions = zeros(4,1);
lowestdepartBound = maxTimeGuess;
lowestapproachBound = maxTimeGuess;
lowestcoastFraction = maxTimeGuess;
lowestYdata = zeros(5,1);
lowestTdata = zeros(5,1);

if direction == Direction.FORWARD
    lowestEscape = NaN;
else
    lowestEscape = maxTimeGuess;
end
%The following for loop will find the optimal MinTimeGuess for the
%current time function so that the optimization algorithm will run more
%quickly


MaxIterations = Constants.FMIN_MAXNUMITER;

refinedTimeStep = largeForwardStep;

foundMinTimeGuess = false;
while ~foundMinTimeGuess %This is simply a check to ensure that we will get a result, do not expect this to run more than once
    %Loop forward at a large interval until the program finds an
    %interesting value of fval
    for timeGuess = lowerTimeLimit : largeForwardStep : upperTimeLimit
        %Randomize coefficients
        cConverge = (pi/2)*(-1 + 2*rand(1,2*order+3));
        if optimize(LegNumber,1)
            departBound(LegNumber,1) = rand(1) * .3 + .051;
            approachBound(LegNumber,1) = rand(1) * .14 + .051;
            coastFraction(LegNumber,1) = rand(1) * .45 + .51;
        end
        
        firstInterestingTime = timeGuess;
        
        [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isInteresting,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
            timeObjectiveFunction(FS,cConverge,order,direction,timeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct, throttleFlag);
        
        %     Store lowest values if converge and if it is lowest trip time
        
        [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
            StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
        
        
        %If an interesting value of fval is found, loop backward from the
        %current timeGuess at a smaller interval until the program finds
        %an uninteresting value of fval
        if (isInteresting == 1)
            fprintf('FOUND FIRST INTERESTING POINT!\n\n')
            for innerTimeGuess = firstInterestingTime + smallerBackStep : smallerBackStep : innerLowerTimeLimit
                
                [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isInteresting,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
                    timeObjectiveFunction(FS,cConverge,order,direction,innerTimeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct, throttleFlag);
                
                %     Store lowest values if converge and if it is lowest trip time
                
                [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
                    StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
                
                
                firstUninterestingTime = innerTimeGuess;
                
                %Once an uninteresting value of fval is found, continue to
                %loop backward until another uninteresting value of fval is
                %found
                if (isInteresting == 0)
                    fprintf('FOUND FIRST UNINTERESTING POINT. STEPPING BACK\n\n')
                    for innerInnerTimeGuess = firstUninterestingTime + smallerBackStep : smallerBackStep : innerLowerTimeLimit
                        
                        [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isInteresting,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
                            timeObjectiveFunction(FS,cConverge,order,direction,innerInnerTimeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct,throttleFlag);
                        
                        %     Store lowest values if converge and if it is lowest trip time
                        
                        [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
                            StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
                        
                        %Once the second uninteresting value of fval is
                        %found, MinTimeGuess is set to the current value of
                        %timeGuess and the entire MinTimeGuess loop is
                        %exited
                        if (isInteresting == 0)
                            fprintf('FOUND MINIMUM TIME GUESS!!\nNOW SEARCHING FOR LOWEST CONVERGENCE\n\n')
                            foundMinTimeGuess = 1;
                            MinTimeGuess = innerInnerTimeGuess;
                            break;
                            
                        end
                    end
                end
                if (foundMinTimeGuess)
                    break;
                end
            end
        end
        if (foundMinTimeGuess)
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following code iterates through time guesses ranging from min to a max.
%This is used instead of the outer fminsearch minimizing the trip time. It
%splits up the process into three parts: the first iterates upward from the
%min time guess and exits when convergence begins to occur. The second loop
%then iterates downward at a smaller time step until convergence
%ceases to occur, at which point it will exit and move on to the third. The
%third loop uses a binary search method to find the minimum trip time which
%must be found between the last converging time and the first
%non-converging time of the previous loop. This whole process requires that
%the min time guess is below the lowest possible trip time

%Initializing the time guess to the min
timeGuess = MinTimeGuess + timeStep;

%First loop
while timeGuess < maxTimeGuess
    %Randomize coefficients
    cConverge = (pi/2)*(-1 + 2*rand(1,2*order+3));
    if optimize(LegNumber,1)
        departBound(LegNumber,1) = rand(1) * .3 + .051;
        approachBound(LegNumber,1) = rand(1) * .14 + .051;
        coastFraction(LegNumber,1) = rand(1) * .45 + .51;
    end
    [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,~,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
        timeObjectiveFunction(FS,cConverge,order,direction,timeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct,throttleFlag);
    
    %     Store lowest values if converge and if it is lowest trip time
    
    [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
        StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
    
    %Exitting when convergence occurs
    if converge
        fprintf('FOUND FIRST CONVERGENCE. STEPPING BACK\n\n')
        break;
    end
    
    %Decreasing by time step specified above
    timeGuess = timeGuess + timeStep;
end


%Approach minimum time guess as closely as possible
multlimit = 3;
for stepMult = 1:1:multlimit
    
    fprintf('STEP MULTIPLIER: %g of %g\n\n',stepMult,multlimit)
    
    %Decrease time guess
    timeGuess = lastConverge - timeStep/(2^stepMult);
    
    while timeGuess > 0
        
        [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,~,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
            timeObjectiveFunction(FS,cConverge,order,direction,timeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct,throttleFlag);
        
        %     Store lowest values if converge and if it is lowest trip time
        
        [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
            StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
        
        %checking to see if this is the flawed part of the code
        
        %         Exiting when convergence ceases to occur twice in a row
        if ~converge
            fprintf('DECREASING TIME STEP\n\n')
            
            refinedTimeStep = timeStep/(2^stepMult);
            %Decreasing by a much smaller time step
            timeGuess = timeGuess - refinedTimeStep;
            
            [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,~,departBound,approachBound,coastFraction,lastYdata,lastTdata] = ...
                timeObjectiveFunction(FS,cConverge,order,direction,timeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct,throttleFlag);
            
            %     Store lowest values if converge and if it is lowest trip time
            
            [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
                StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata);
            
            %Exit if does not converge
            if ~converge
                break
            end
            
        end
        refinedTimeStep = timeStep/(2^stepMult);
        %Decreasing time guess
        timeGuess = timeGuess - refinedTimeStep;
    end
    
end

if direction == Direction.RETURN
    lowestReturnConverge = lowestConverge;
    lowestForwardConverge = nan;
else
    lowestForwardConverge = lowestConverge;
    lowestReturnConverge = nan;
end

%Call IntegrateForPlot to integrate lowest values for plotting
[FuelMass,Ydata,Tdata,Adata,departData,coastData,thrustData,approachData,MdotData,TmdotData] = IntegrateForPlot(COAST,direction,LegNumber,efficiency,mdot,numberOfEngines,inputPower,lowestEscape,lowestfinalMass,numIntervals,interpolateStruct,throttleFlag,lowestYdata,lowestTdata);

finalTime = lowestConverge;
escapeTime = lowestEscape;
phiCoeff = lowestcConverge;
fval = lowestfval;

asteroidConditions = lowestasteroidConditions;

%Plot spacecraft motion variables with respect to time
[plotFileName,rError,thetaError,uError,hError] = ...
    PlotAndWriteToFile(order,direction,phiCoeff,fval,escapeTime,finalTime,timeWait,uExhaust,asteroidCloseApp,asteroidConditions,astDes,powerLevel,thrusterName,Tdata,Ydata,thrusterStructs,Adata,departData,coastData,thrustData,approachData,COAST,LegNumber,AsteroidSemimajor,AsteroidEccentricity,MdotData,TmdotData);

%Ends timer for function
timeToExecute = toc;
end
