%Main function; iterates through chosen thrusters and powers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function contains variables for minumum and maximum values that are
%passed to the MinimizeTripTime function. It also contains variables for
%incrementing those values. The function will iterate MinimizeTripTime over
%the given range of parameters with the specified increment size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Main function of program
function MultipleRuns()

clear all
format compact
clear classes

%Prints the Console to a text file (Edit out if not using overnight as well
%as diary off at end of MultipleRuns)
ConsoleOutputName = [freename('./','ConsoleOutput') '.txt'];
diary (ConsoleOutputName)

%Let user decide whether to read input from text file or manual input trip
%info

%Seeding the random number generator so that the same
%results are not produced each time the rand function is
%called
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

%Let user decide whether to read data from the Excel file or input manually
dataEntry = input('Would you like to manually enter specifications (0) or read from Input.xlsx file(1)? ');

%Initialize profiler decision
pfl = 0;

if dataEntry ~= 0
    %The following starts the profiler if desired from the excel file.
    pfl = xlsread('Input.xlsx','Base Inputs','c1:c1');
    if pfl ~= 0
        profile on
        fprintf('Profiler is On\n')
    end
end

%%%%%%%%%%%%%%%%Hard Coded Info%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set order for Fourier Series evaluation
minOrder = 7;
maxOrder = 7;

% [flagStore] = ChooseMinimize();
flagStore = true;

%Mass of the structure of the spaceship
%Includes mass of fuel tank
massStruct = 0;

%Mass of the sample of the asteroid being brought back. Only used in the
%return journey.
massSample = 10;

%Allow user to choose starting altitude of mission
%selectForwardAltitude = input('Input starting altitude (m) measured from the surface of earth for mission: ');
selectForwardAltitude = 500000; %This is hard coded right now

%Allow user to choose ending altitude of mission
%selectReturnAltitude = input('Input ending altitude (m) measured from the surface of earth for mission: ');
selectReturnAltitude = 120000; %This is hard coded right now

%Magnitude of velocity at returnRadius (m/s)
%Made up of radial velocity and angular velocity components:
%sqrt(u^2+(h/r)^2)
%entryVelocity = input('Input entry velocity to earth (m/s) at ending altitude: ');
entryVelocity = 12000; %This is hard coded right now

%Entry angle at returnRadius (deg)
%This angle is the inverse tangent of radial velocity over angular
%velocity: atan(u*r/h)
%entryAngle = input('Input entry angle to earth (m/s) at ending altitude: ');
entryAngle = 7; %This is hard coded right now

%Power input information

%SOLAR POWER VS. NUCLEAR INPUT OPTION
isSolar = 1; %0 for nuclear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dataEntry == 0
    %Number of trials per thruster
    exitFlag = false;
    while ~exitFlag
        exitFlag = true;
        numTrials = input('Enter number of trials per thruster: ');
        if numTrials <= 0 || mod(numTrials,1) ~= 0
            fprintf('\nFatal Input. Please enter an integer greater than 0.\n\n')
            exitFlag = false;
        end
    end
    
    %Allow user to choose coasting profile
    [COAST,departBound,approachBound,numIntervals,coastFraction,optimize] = ChooseCoast(dataEntry);
    
    %Allow user to choose alphaP values
    [minAlphaP,maxAlphaP] = ChooseAlphaP();
    
    %Allow user to choose payload of mission
    massPayload = ChoosePayload();
    
    %Allow user to choose asteroid
    [astDes] = ChooseAsteroid(dataEntry);
    
    %Get orbital period of chosen asteroid
    [~,~,AsteroidEccentricity,AsteroidSemimajor,AsteroidPeriod,zAxisFraction] = ParseFile('AsteroidPosition',astDes);
    
    % %Display Oribital Period for user reference
    % fprintf('The orbital period of this asteroid is: %g\n', AsteroidPeriod)
    
    %Allow user to choose length of stay on asteroid
    [N] = input('Stay time on the asteroid (in half-orbital periods): ');
    
    %Allow user to choose between fully converged or quick run.  This is used
    %for testing purposes.  If you would like to make sure the code runs
    %without errors, use a fast trial. If you would like the best data use
    %fully converged.
    [FS] = input('Fast trial(0) or Fully Converged(1)? ');
    
    %Gives the user the choice to use the current know data that is
    %provided by the manufacturer or to have GetInterpolData extrapolate
    %the data.
    chooseExtrap = input('Use tested Pin and Uex data(0) or Extrapolate data(1)? ');
    if chooseExtrap == 0
        extrapPercent = 0;
    elseif chooseExtrap == 1 
        extrapPercent = input('Enter the percent of extrapolation (in decimal form) :');
    end
    
    
     %Give user option to specify preference of higher thrust (0) or efficiency (1) 
    throttleFlag = input('Preference of higher thrust (0) or efficiency (1): ');
    
    
    
    %Calculate time wait based on user input
    if N == 0
        timeWait = 0;
    elseif mod(N,2) ~= 0
        timeWait = zAxisFraction * AsteroidPeriod + ((N-1)/2) * AsteroidPeriod;
    elseif mod(N,2) == 0
        timeWait = (1/2)* N * AsteroidPeriod;
    end
    
    fprintf('Time Wait: %g\n',timeWait)
    
    %Allow user to choose any combination of thrusters
    [numChoiceThrusters,thrusterStructs,nameofThruster] = ChooseThrusters(dataEntry, throttleFlag);
    
    %Initialize the lowest values based on the number of thrusters that the
    %user inputted
    lowTrial = zeros(1,numChoiceThrusters);
    lowTime = zeros(1,numChoiceThrusters);
    lowestforwardEscapeTime = zeros(1,numChoiceThrusters);
    lowestforwardTripTime = zeros(1,numChoiceThrusters);
    lowestWaitTime = zeros(1,numChoiceThrusters);
    lowestreturnTripTime = zeros(1,numChoiceThrusters);
    lowestreturnEscapeTime = zeros(1,numChoiceThrusters);
    lowestMdot = zeros(1,numChoiceThrusters);
    lowestPPMass = zeros(1,numChoiceThrusters);
    lowestFuelMass = zeros(1,numChoiceThrusters);
    lowestWetMass = zeros(1,numChoiceThrusters);
    lowestTotalTimeError = zeros(1,numChoiceThrusters);
    
    %Allow user to choose either one or all power levels
    [powerStart,powerEnd] = ChoosePowerLevels();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %Number of trials per thruster
    exitFlag = false;
    while ~exitFlag
        exitFlag = true;
        numTrials = xlsread('Input.xlsx','Base Inputs','b1:b1');
        if numTrials <= 0 || mod(numTrials,1) ~= 0
            fprintf('\nFatal Input. Please enter an integer greater than 0.\n\n')
            exitFlag = false;
        end
    end
    fprintf('Number of Trials per thruster: %g\n',numTrials)
    
    %Sets coasting info from data in excel file
    [COAST,departBound,approachBound,numIntervals,coastFraction,optimize] = ChooseCoast(dataEntry);
    fprintf('Coast during forward trip? : %g\n',COAST(2,1));
    if COAST(2,1)
        fprintf('Coasting intervals during forward trip = %g\n',numIntervals(2,1));
    end
    if ~isnan(departBound(2,1)) && ~isnan(approachBound(2,1)) && ~isnan(coastFraction(2,1))
        fprintf('Coasting fraction per interval during forward trip = %g\n',coastFraction(2,1));
        fprintf('Depart Bound during forward trip = %g\n',departBound(2,1));
        fprintf('Approach Bound during forward trip = %g\n',approachBound(2,1));
    end
    fprintf('Coast during return trip? : %g\n',COAST(1,1));
    if COAST(1,1)
        fprintf('Coasting intervals during return trip = %g\n',numIntervals(1,1));
    end
    if ~isnan(departBound(1,1)) && ~isnan(approachBound(1,1)) && ~isnan(coastFraction(1,1))
        fprintf('Coasting fraction per interval during return trip = %g\n',coastFraction(1,1));
        fprintf('Depart Bound during return trip = %g\n',departBound(1,1));
        fprintf('Approach Bound during return trip = %g\n',approachBound(1,1));
    end
    
    %Set alphaP values from data in excel file
    minAlphaP = xlsread('Input.xlsx','Base Inputs','b2:b2');
    maxAlphaP = xlsread('Input.xlsx','Base Inputs','b3:b3');
    fprintf('Minimum AlphaP : %g\n',minAlphaP)
    fprintf('Maximum AlphaP : %g\n',maxAlphaP)
    
    %Set payload values from data in excel file
    massPayload = xlsread('Input.xlsx','Base Inputs','b4:b4');
    fprintf('Mass Payload : %g',massPayload)
    %Set asteroid data based on excel file
    [astDes,astChoice] = ChooseAsteroid(dataEntry);
    fprintf('Asteroid Designation : %g %s\n',astChoice,astDes)
    
    %Get orbital period of chosen asteroid
    [~,~,AsteroidEccentricity,AsteroidSemimajor,AsteroidPeriod,zAxisFraction] = ParseFile('AsteroidPosition',astDes);
    
    % %Display Oribital Period for user reference
    % fprintf('The orbital period of this asteroid is: %g\n', AsteroidPeriod)
    
    %From excel file, input choice between fully converged or quick run.  This is used
    %for testing purposes.  If you would like to make sure the code runs
    %without errors, use a fast trial. If you would like the best data use
    %fully converged.
    FS = xlsread('Input.xlsx','Base Inputs','b7:b7');
    fprintf('Fast Convergence (0) or Fully Converged(1) : %g\n',FS)
    
    %This shows the user whether or not the function will extrapolate the
    %data that is given by the manufacturer.
    chooseExtrap = xlsread('Input.xlsx','Operation Guidelines','b17:b17');
    fprintf('Use tested Pin and Uex data(0) or Extrapolate data(1) : %g\n',chooseExtrap)
    
    extrapPercent = xlsread('Input.xlsx','Operation Guidelines','b18:b18');
    fprintf('Percent of extrapolation: %g\n', extrapPercent)
    
      %read from excel preference of higher thrust or higher efficiency
    throttleFlag = xlsread('Input.xlsx','Operation Guidelines','b19:b19');
    fprintf('Preference of higher thrust (0) or efficiency (1) : %g\n',throttleFlag);
    
    
    %read in stay time from excel file
    N = xlsread('Input.xlsx','Base Inputs','b6:b6');
    fprintf('Stay time on asteroid (in half-orbital periods) : %g\n',N)
    
    %Calculate time wait based on excel input
    if N == 0
        timeWait = 0;
    elseif mod(N,2) ~= 0
        timeWait = zAxisFraction * AsteroidPeriod + ((N-1)/2) * AsteroidPeriod;
    elseif mod(N,2) == 0
        timeWait = (1/2)* N * AsteroidPeriod;
    end
    
    fprintf('Time Wait: %g\n',timeWait)
    
    %Read in combination of thrusters
    [numChoiceThrusters,thrusterStructs,nameofThruster] = ChooseThrusters(dataEntry, throttleFlag);
    fprintf('Number of thrusters to be evaluated : %g\n',numChoiceThrusters)
    fprintf('Thrusters being evaluated: \n')
    for i=1:1:numChoiceThrusters
        fprintf('%s \n',nameofThruster{i})
    end
    %Initialize the lowest values based on the number of thrusters that the
    %user inputted
    lowTrial = zeros(1,numChoiceThrusters);
    lowTime = zeros(1,numChoiceThrusters);
    lowestforwardEscapeTime = zeros(1,numChoiceThrusters);
    lowestforwardTripTime = zeros(1,numChoiceThrusters);
    lowestWaitTime = zeros(1,numChoiceThrusters);
    lowestreturnTripTime = zeros(1,numChoiceThrusters);
    lowestreturnEscapeTime = zeros(1,numChoiceThrusters);
    lowestMdot = zeros(1,numChoiceThrusters);
    lowestPPMass = zeros(1,numChoiceThrusters);
    lowestFuelMass = zeros(1,numChoiceThrusters);
    lowestWetMass = zeros(1,numChoiceThrusters);
    lowestTotalTimeError = zeros(1,numChoiceThrusters);
    
    %Read power levels from the excel file
    powerStart = xlsread('Input.xlsx','Base Inputs','b8:b8');
    powerEnd = xlsread('Input.xlsx','Base Inputs','b9:b9');
    fprintf('Lower engine multiplier: %g\n',powerStart)
    fprintf('Upper engine multiplier: %g\n',powerEnd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%Open file for writing
Data = fopen('MultipleRunsData.txt','w');
fprintf(Data, 'Direction\tThruster\tAsteroid\tTotal Power\tSolar Power\tWait Time\tTotal Time(yr)\tLeg Time\tEscape Time\tComputation Time(s)\tEarth Escape Angle(rad)\tFinal Function Value\tNumber of Engines\tAlphaP\tAlphaT\tInput Power\tUex\tEfficiency\tR error\tTheta error\tU error\tH error\tPlot File Name\tMdot\tPayloadMass\tStruct Mass\tPPlant Mass\tFuel Mass\tThruster Mass\tSample Mass\tWet Mass\r\n');
%Print data headers to the file.
%Note the double tab between the thruster angle coefficients and the rest
%of the data. This is intentional, for future importing into Excel.
%Close file after writing
fclose(Data);

%Iterate over polynomial-fit order
for order = minOrder:1:maxOrder
    
    if maxAlphaP > minAlphaP
        dAp = (maxAlphaP - minAlphaP);
    else
        dAp = 16;
    end
    
    %Iterate through alpha P values
    for alphaP = minAlphaP:dAp:maxAlphaP
        
        %Iterate through the array of user selected thrusters
        for n = 1:1:numChoiceThrusters
            [interpolateStruct] = GetInterpolData(thrusterStructs(n).InterpolDataFileName,thrusterStructs(n).ThrusterIDNumber,chooseExtrap,extrapPercent);
            %set max efficiency, and Uex and Mdot for final output value
            if throttleFlag == 0
                efficiency = interpolateStruct.etaConstMDOT;
                Uex = interpolateStruct.uExConstMDOT;
                mDot = interpolateStruct.mDotConstMDOT;
                maxEfficiency = (efficiency(1,1));
                maxUex = Uex(1,1);
                maxMdot = mDot(1,1);
                
            else
                efficiency = interpolateStruct.etaConstUEX;
                Uex = interpolateStruct.uExConstUEX;
                mDot = interpolateStruct.mDotConstUEX;
                maxEfficiency = (efficiency(1,1));
                maxUex = Uex(1,1);
                maxMdot = mDot(1,1);
                
            end
            
            
            
            
            %Iterate through the array of user selected power levels
            for e = powerStart:1:powerEnd
                
                %Iterate through several runs of each case
                for numPoints = 1:1:numTrials
                    %Return journey is evaluated first in order to
                    %determine first the minimum amount of fuel needed to
                    %return. This fuel is added onto the forward trip
                    
                    %Set to return journey.
                    direction = Direction.RETURN;
                    %Initialize the total trip time
                    timeRoundTrip = 0;
                    %Initialize the total mass of fuel
                    totalMassFuel = 0;
                    %Mass of fuel (kg)
                    massFuel = 0;
                    %Time of last convergence
                    lastConverge = 0;
                    
                    %Total power for a thruster is the number of engines used
                    %multiplied by its input power.
                    totalPower = (thrusterStructs(n).NumEnginesMultiplier * e) * (thrusterStructs(n).InputPower + (thrusterStructs(n).InputPower * extrapPercent));
                    
                    %This loop runs through the return journey and then the
                    %forward journey and calculates the total round trip
                    %time.
                    for LegNumber = 1:1:2
                        %1 is return, 2 is forward
                        %Shows user whether MATLAB is currently calculating
                        %the forward journey or the return journey
                        if LegNumber==1
                            fprintf('____________________________________________________________________________________\n')
                            fprintf('                            %s Total Power: %g                                      \n',thrusterStructs(n).ThrusterName,totalPower)
                            fprintf('                                    Trial: %g                                       \n',numPoints)
                            fprintf('                                                                                    \n')
                            fprintf('                                   CALCULATING:                                     \n')
                            fprintf('                                    ~RETURN~                                        \n')
                            fprintf('                                    ~JOURNEY~                                       \n')
                            fprintf('____________________________________________________________________________________\n')
                        elseif LegNumber==2
                            fprintf('____________________________________________________________________________________\n')
                            fprintf('                            %s Total Power: %g                                      \n',thrusterStructs(n).ThrusterName,totalPower)
                            fprintf('                                    Trial: %g                                       \n',numPoints)
                            fprintf('                                                                                    \n')
                            fprintf('                                   CALCULATING:                                     \n')
                            fprintf('                                    ~FORWARD~                                       \n')
                            fprintf('                                    ~JOURNEY~                                       \n')
                            fprintf('____________________________________________________________________________________\n')
                        end
                        %Base initial time guess off of the previous trial
                        if numPoints == 1;
                            maxTimeLimit = nan;
                        else
                            if direction == Direction.RETURN
                                maxTimeLimit = lowestReturnConverge;
                            else
                                maxTimeLimit = lowestForwardConverge;
                            end
                        end
                        
                        %Call MinimizeTripTime
                        %MinimizeTripTimeFast only does forward time step
                        %until convergence.
                        %MinimizeTripTimeFull looks for fully converged
                        %case
                        %FS is choice of fast convergence or full
                        %convergence
                        
                        if FS == 0
                            
                            if direction == Direction.RETURN
                                [finalTime,escapeTime,runTime,phiCoeff,plotFileName,rError,thetaError,uError,hError,fVal,lastConverge,massPowerPlant,massThruster,massStruct,mdot,lowestEscape,lowestConverge,lowestReturnConverge,~,returnTimeError,massFuel,departBound,approachBound,coastFraction] = ...
                                    MinimizeTripTimeFast(FS,order,direction,massFuel,(thrusterStructs(n).NumEnginesMultiplier * e),alphaP,thrusterStructs(n).AlphaT,thrusterStructs(n).InputPower,thrusterStructs(n).Uex,thrusterStructs(n).Efficiency,massPayload,massStruct,massSample,lastConverge,selectForwardAltitude,selectReturnAltitude,entryVelocity,entryAngle,astDes,timeWait,isSolar,totalPower,thrusterStructs(n).ThrusterName,numPoints,thrusterStructs(n),maxTimeLimit,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,flagStore,optimize,interpolateStruct, throttleFlag);
                            else
                                [finalTime,~,runTime,phiCoeff,plotFileName,rError,thetaError,uError,hError,fVal,lastConverge,massPowerPlant,massThruster,massStruct,mdot,~,lowestConverge,~,lowestForwardConverge,forwardTimeError,massFuel,departBound,approachBound,coastFraction] = ...
                                    MinimizeTripTimeFast(FS,order,direction,massFuel,(thrusterStructs(n).NumEnginesMultiplier * e),alphaP,thrusterStructs(n).AlphaT,thrusterStructs(n).InputPower,thrusterStructs(n).Uex,thrusterStructs(n).Efficiency,massPayload,massStruct,massSample,lastConverge,selectForwardAltitude,selectReturnAltitude,entryVelocity,entryAngle,astDes,timeWait,isSolar,totalPower,thrusterStructs(n).ThrusterName,numPoints,thrusterStructs(n),maxTimeLimit,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,flagStore,optimize,interpolateStruct, throttleFlag);
                            end
                        else
                            if direction == Direction.RETURN
                                [finalTime,escapeTime,runTime,phiCoeff,plotFileName,rError,thetaError,uError,hError,fVal,lastConverge,massPowerPlant,massThruster,massStruct,mdot,lowestEscape,lowestConverge,lowestReturnConverge,~,returnTimeError,massFuel,departBound,approachBound,coastFraction] = ...
                                    MinimizeTripTimeFull(FS,order,direction,massFuel,(thrusterStructs(n).NumEnginesMultiplier * e),alphaP,thrusterStructs(n).AlphaT,thrusterStructs(n).InputPower,thrusterStructs(n).Uex,thrusterStructs(n).Efficiency,massPayload,massStruct,massSample,lastConverge,selectForwardAltitude,selectReturnAltitude,entryVelocity,entryAngle,astDes,timeWait,isSolar,totalPower,thrusterStructs(n).ThrusterName,numPoints,thrusterStructs(n),maxTimeLimit,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,flagStore,optimize,interpolateStruct, throttleFlag);
                            else
                                [finalTime,~,runTime,phiCoeff,plotFileName,rError,thetaError,uError,hError,fVal,lastConverge,massPowerPlant,massThruster,massStruct,mdot,~,lowestConverge,~,lowestForwardConverge,forwardTimeError,massFuel,departBound,approachBound,coastFraction] = ...
                                    MinimizeTripTimeFull(FS,order,direction,massFuel,(thrusterStructs(n).NumEnginesMultiplier * e),alphaP,thrusterStructs(n).AlphaT,thrusterStructs(n).InputPower,thrusterStructs(n).Uex,thrusterStructs(n).Efficiency,massPayload,massStruct,massSample,lastConverge,selectForwardAltitude,selectReturnAltitude,entryVelocity,entryAngle,astDes,timeWait,isSolar,totalPower,thrusterStructs(n).ThrusterName,numPoints,thrusterStructs(n),maxTimeLimit,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,flagStore,optimize,interpolateStruct, throttleFlag);
                            end
                        end
                        
                        %Update the total mass of fuel
                        totalMassFuel = totalMassFuel + massFuel;
                        %Update the total wet mass of spacecraft
                        wetMass = totalMassFuel + massStruct + massPayload + massPowerPlant + massThruster;
                        %Update the total round trip time
                        timeRoundTrip = timeRoundTrip + lowestConverge;
                        
                        %Save results for each direction
                        if direction==Direction.RETURN
                            returnTripTime=lowestConverge;
                            returnEscapeTime = lowestEscape;
                            returnRunTime = runTime;
                            returnDepartBound = departBound(1,1);
                            returnApproachBound = approachBound(1,1);
                            returnCoastFraction = coastFraction(1,1);
                            returnNumIntervals = numIntervals(1,1);
                        elseif direction==Direction.FORWARD
                            forwardTripTime=lowestConverge;
                            forwardRunTime = runTime;
                            forwardDepartBound = departBound(2,1);
                            forwardApproachBound = approachBound(2,1);
                            forwardCoastFraction = coastFraction(2,1);
                            forwardNumIntervals = numIntervals(2,1);
                        end
                        
                        if direction == Direction.RETURN
                            tripDirection = 'Return';
                        elseif direction == Direction.FORWARD
                            tripDirection = 'Forward';
                        end
                        
                        %Mass of fuel for return trip is now known, so
                        %the forward trip can now be evaluated.
                        %Set for forward journey
                        direction = Direction.FORWARD;
                        
                        %Open file for appending
                        Data = fopen('MultipleRunsData.txt','a');
                        
                        %Outputs important imformation into data file in a tabular format
                        %in the order to be imported into Excel
                        fprintf(Data, '%s\t%s\t%s\t%6.12g\t%i\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%s\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\t%6.12g\r\n',...
                            tripDirection,thrusterStructs(n).ThrusterName,astDes,totalPower,isSolar,timeWait,timeRoundTrip,finalTime,escapeTime,runTime,phiCoeff(end),fVal,(thrusterStructs(n).NumEnginesMultiplier * e),alphaP,thrusterStructs(n).AlphaT,thrusterStructs(n).InputPower,thrusterStructs(n).Uex,thrusterStructs(n).Efficiency(1),rError,thetaError,uError,hError,plotFileName,mdot,massPayload,massStruct,massPowerPlant,massFuel,massThruster,massSample,wetMass);
                        %Close file
                        fclose(Data);
                    end
                    
                    %Sum total mission time, associated errors and total
                    %compuation time.
                    totalMissionTime = forwardTripTime + timeWait + returnTripTime + returnEscapeTime;
                    totalTimeError = sqrt(returnTimeError^2 + forwardTimeError^2);
                    totalRunTime = returnRunTime + forwardRunTime;
                    
                    %Output each leg of the trip and the total mission time
                    fprintf(' ____________________________________________________________________________________\n')
                    fprintf('              %s - %s - %g kW                                                        \n',astDes,thrusterStructs(n).ThrusterName,totalPower)
                    fprintf('                                    Trial: %g                                        \n', numPoints)
                    fprintf('                                                                                     \n')
                    fprintf('                             FORWARD TIME: %g                                        \n',forwardTripTime)
                    fprintf('                                WAIT TIME: %g (%g Half Orbital Periods)              \n',timeWait, N)
                    fprintf('                              RETURN TIME: %g                                        \n',returnTripTime)
                    fprintf('                       RETURN ESCAPE TIME: %g                                        \n',returnEscapeTime)
                    fprintf('                                          **************                             \n')
                    fprintf('                          TOTAL TRIP TIME: %g     +/-  %g                            \n',totalMissionTime,totalTimeError)
                    fprintf('_____________________________________________________________________________________\n')
                    
                    if ~COAST(1,1)
                        returnDepartBound = NaN;
                        returnApproachBound = NaN;
                        returnCoastFraction = NaN;
                        returnNumIntervals = NaN;
                    end
                    
                    if ~COAST(2,1)
                        forwardDepartBound = NaN;
                        forwardApproachBound = NaN;
                        forwardCoastFraction = NaN;
                        forwardNumIntervals = NaN;
                    end
                    %Assume that the first run was the lowest total mission
                    %time and then assign the corresponding times to the
                    %lowest times
                    if numPoints == 1
                        lowTrial(1,n) = numPoints;
                        lowTime(1,n) = totalMissionTime;
                        lowestforwardTripTime(1,n) = forwardTripTime;
                        lowestWaitTime(1,n) = timeWait;
                        lowestreturnTripTime(1,n) = returnTripTime;
                        lowestreturnEscapeTime(1,n) = returnEscapeTime;
                        lowestMdot(1,n) = mdot;
                        lowestPPMass(1,n) = massPowerPlant;
                        lowestFuelMass(1,n) = totalMassFuel;
                        lowestWetMass(1,n) = totalMassFuel(1,n) + lowestPPMass(1,n) + massStruct + massPayload + massThruster;
                        lowestTotalTimeError(1,n) = totalTimeError;
                        lowestReturnCoastFraction(1,n) = returnCoastFraction;
                        lowestReturnApproachBound(1,n) = returnApproachBound;
                        lowestReturnDepartBound(1,n) = returnDepartBound;
                        lowestForwardCoastFraction(1,n) = forwardCoastFraction;
                        lowestForwardApproachBound(1,n) = forwardApproachBound;
                        lowestForwardDepartBound(1,n) = forwardDepartBound;
                        %Compare lowest time stored to the total mission time
                        %produced during this loop and then override values if
                        %it is smaller
                    elseif totalMissionTime < lowTime(1,n)
                        lowTrial(1,n) = numPoints;
                        lowTime(1,n) = totalMissionTime;
                        lowestforwardTripTime(1,n) = forwardTripTime;
                        lowestWaitTime(1,n) = timeWait;
                        lowestreturnTripTime(1,n) = returnTripTime;
                        lowestreturnEscapeTime(1,n) = returnEscapeTime;
                        lowestMdot(1,n) = mdot;
                        lowestPPMass(1,n) = massPowerPlant;
                        lowestFuelMass(1,n) = totalMassFuel;
                        lowestWetMass(1,n) = totalMassFuel(1,n) + lowestPPMass(1,n) + massStruct + massPayload + massThruster;
                        lowestTotalTimeError(1,n) = totalTimeError;
                        lowestReturnCoastFraction(1,n) = returnCoastFraction;
                        lowestReturnApproachBound(1,n) = returnApproachBound;
                        lowestReturnDepartBound(1,n) = returnDepartBound;
                        lowestForwardCoastFraction(1,n) = forwardCoastFraction;
                        lowestForwardApproachBound(1,n) = forwardApproachBound;
                        lowestForwardDepartBound(1,n) = forwardDepartBound;
                    end
                    
                    
                    gamma = mod(phiCoeff(end),2*pi) * (180/(2 * pi^2));
                    theta = mod(phiCoeff(end-1),2*pi) * (180/(2 * pi^2));
                    
                    %Create seperate excel files for each trial
                    TrialFileName = freename('./',thrusterStructs(n).ThrusterName,3);
                    %Format relevent data
                    TrialSummaryData = {'','','',astDes,'','','','';
                        '','','',thrusterStructs(n).ThrusterName,'','','','';
                        '','','Trial:',numPoints,'','','','';
                        '','','','','','','','';
                        'Trial Summary:','','','','','','','';
                        '','Theta:',theta,'','Gamma:',gamma,'','';
                        'Total Trip Time:',totalMissionTime,'+/-',totalTimeError,'','','','';
                        'Time by Leg','','','','','','','';
                        '','Forward Helio','Wait','Return Helio','Return Escape','','','';
                        '',forwardTripTime,timeWait,returnTripTime,returnEscapeTime,'','','';
                        'Coasting Details','','','','','','','';
                        'Forward:','','','','','','','';
                        '','Depart Fraction','Coasting Segment','Approach Fraction','','','','';
                        '',forwardDepartBound,1 - forwardDepartBound - forwardApproachBound,forwardApproachBound,'','','','';
                        '','Num Intervals','Coast Fraction','Thrust Fraction','','','','';
                        '',forwardNumIntervals,forwardCoastFraction,1 - forwardCoastFraction,'','','','';
                        'Return:','','','','','','','';
                        '','Depart Fraction','Coasting Segment','Approach Fraction','','','','';
                        '',returnDepartBound,1 - returnDepartBound - returnApproachBound,returnApproachBound,'','','','';
                        '','Num Intervals','Coast Fraction','Thrust Fraction','','','','';
                        '',returnNumIntervals,returnCoastFraction,1 - returnCoastFraction,'','','','';
                        'Engine Details','','','','','','','';
                        'Name','Power','Engine Multiplier','Total Power','Efficiency','Alpha P','Alpha T','Mdot';
                        thrusterStructs(n).ThrusterName,(thrusterStructs(n).InputPower + (thrusterStructs(n).InputPower * extrapPercent)),powerStart,totalPower,maxEfficiency,alphaP,thrusterStructs(n).AlphaT,maxMdot;
                        'Payload Mass','Structure Mass','PPlant Mass','Fuel Mass','Thruster Mass','Sample Mass','Wet Mass','';
                        massPayload,massStruct,massPowerPlant,totalMassFuel,massThruster,massSample,wetMass,'';
                        'Computation Time:',totalRunTime,'','','uExhaust:',maxUex,'','';
                        'Forward:',forwardRunTime,'','Return:',returnRunTime,'','','';
                        'Fast/Slow Converge','','','ThrottleFlag','','','Extrap Percent','';
                        FS,'','',throttleFlag,'','',extrapPercent,'';};
                    %Export to excel
                    xlswrite(TrialFileName,TrialSummaryData,numPoints)
                end
                
                %Check for ease of comparison while running
                fprintf('LOW TIME: %g\n', lowTime(1,n))
                
                %Format relevent data of best trial per engine multipler
                MissionSummaryData = {'','','',astDes,'','','','';
                    '','','',thrusterStructs(n).ThrusterName,'','','','';
                    '','','Trial:',lowTrial(1,n),'','','','';
                    '','','','','','','','';
                    'Important Info:','Fuel Mass','Wet Mass','Trip Time','','','','';
                    '',lowestFuelMass(1,n),lowestWetMass(1,n),lowTime(1,n),'','','','';
                    '','','','','','','','';
                    'Mission Summary:','','','','','','','';
                    'Total Trip Time:',lowTime(1,n),'+/-',lowestTotalTimeError(1,n),'','','','';
                    'Time by Leg','','','','','','','';
                    '','','Forward Helio','Wait','Return Helio','Return Escape','','';
                    '','',lowestforwardTripTime(1,n),timeWait,lowestreturnTripTime(1,n),lowestreturnEscapeTime(1,n),'','';
                    '','','','','','','','';
                    'Coasting Details','','','','','','','';
                    '','','','','','','','';
                    'Forward:','','','','','','','';
                    '','Depart Fraction','Coasting Segment','Approach Fraction','','','','';
                    '',lowestForwardDepartBound,1 - lowestForwardDepartBound - lowestForwardApproachBound,lowestForwardApproachBound,'','','','';
                    '','Num Intervals','Coast Fraction','Thrust Fraction','','','','';
                    '',forwardNumIntervals,lowestForwardCoastFraction,1 - lowestForwardCoastFraction,'','','','';
                    '','','','','','','','';
                    'Return:','','','','','','','';
                    '','Depart Fraction','Coasting Segment','Approach Fraction','','','','';
                    '',lowestReturnDepartBound,1 - lowestReturnDepartBound - lowestReturnApproachBound,lowestReturnApproachBound,'','','','';
                    '','Num Intervals','Coast Fraction','Thrust Fraction','','','','';
                    '',returnNumIntervals,lowestReturnCoastFraction,1 - lowestReturnCoastFraction,'','','','';
                    '','','','','','','','';
                    'Engine Details','','','','','','','';
                    '','','','','','','','';
                    'Name','Power','Engine Multiplier','Total Power','Efficiency','Alpha P','Alpha T','Mdot';
                    thrusterStructs(n).ThrusterName,(thrusterStructs(n).InputPower + (thrusterStructs(n).InputPower * extrapPercent)),powerStart,totalPower,maxEfficiency,alphaP,thrusterStructs(n).AlphaT,maxMdot;
                    'Uexhaust:',maxUex,'','','','','','';
                    'Payload Mass','Structure Mass','PPlant Mass','Fuel Mass','Thruster Mass','Sample Mass','Wet Mass','';
                    massPayload,massStruct,lowestPPMass(1,n),lowestFuelMass(1,n),massThruster,massSample,lowestWetMass(1,n),'';
                    'Fast/Slow Converge','','','ThrottleFlag','','','Extrap Percent','';
                    FS,'','',throttleFlag,'','',extrapPercent,'';};
                %Create seperate excel files for each engine mult & alpha P
                aP = int2str(alphaP);
                fn = strcat(aP,thrusterStructs(n).ThrusterName);
                FileName = freename('./',fn,2);
                %Export to Excel
                xlswrite(FileName,MissionSummaryData,e)
            end
        end
    end
end
diary off;

if dataEntry ~= 0 && pfl ~= 0
    
    %The following finishes the profile and saves it.  To view, click on the
    %Profiler0XX.mat file in matlab and enter the command:
    %profview(0,p)
    %This will bring up the profile as if a run and time had been done.
    p=profile('info');
    ProfileName = [freename('./','Profiler') '.mat'];
    save(ProfileName,'p')
    profile off;
    clear p;
    profile clear;
end
end