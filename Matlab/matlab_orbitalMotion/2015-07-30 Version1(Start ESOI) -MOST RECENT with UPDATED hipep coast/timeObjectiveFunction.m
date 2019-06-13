%TIMEOBJECTIVEFUNCTION Determines if time guess is adequate based on F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function receives the trip time guess variable and uses fminsearch to
%determine whether this is an adequate guess based on the value of F.
%fminsearch will optimize the trip parameters, attempting to land on the
%asteroid. Time minimization occurs by the progressively lower time guesses
%in MinimizeTripTime. Normally, this task would be handled by a
%constrained optimization algorithm but fminsearch is not designed to
%handle constraints. Our solution was to artificially create the global
%minimum of our function at that point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function [cConverge,converge,lastConverge,lastfval,lastEscape,mdot,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isInteresting,departBound,approachBound,coastFraction,y,t] = timeObjectiveFunction(FS,cConverge,order,direction,timeGuess,timeWait,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample,forwardRadius,returnRadius,entryVelocity,entryAngle,lastConverge,lastfval,lastEscape,earthCloseApp,asteroidCloseApp,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,LegNumber,MaxIterations,optimize,interpolateStruct,throttleFlag)


NUMBER_OF_INDICIES = Constants.NUMBER_OF_INDICIES; %initialize constants

%Displays the current time guess on the screen
fprintf('TimeGuess  Value: %g\n',timeGuess)

%Set the time limit for ode45 integration in EscapeEarth
maxEscapeTime = 10;

%Given the ship parameters and user specified initial conditions (all
%global variables), EscapeEarth calculates the time required to leave earth
%sphere of influence. The time passed to it must be greater than the time
%required to escape or the value returned will be incorrect.
if direction == Direction.RETURN
    %Calculates engine parameters given the time guess
    [mdot,accelFinal,finalMass] = CalculateEngineParameters(direction,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample);
    [escapeTime,escapeVelocity,escapeEarthAngle,~,~] = EscapeEarth(order,direction,maxEscapeTime,uExhaust,accelFinal,forwardRadius,returnRadius,entryVelocity,entryAngle,mdot,numberOfEngines);
    
elseif direction == Direction.FORWARD
    [mdot,~,finalMass] = CalculateEngineParameters(direction,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample);
    escapeTime = NaN;
    escapeVelocity = NaN;
    escapeEarthAngle = NaN;
end

%If evaluating return trip, asteroid condiditions must be reevaluated based
%on the expected wait time on the asteroid
if direction == Direction.RETURN
    %Element is 0 to signify asteroid
    orbitalElement = 0;
    
    %Calculates conditions of the asteroid at time of takeoff from the
    %asteroid after waiting for the certain amount of time. The time passed
    %to it is the wait time. If wait is 0, use close approach data.
    if timeWait > 0
        asteroidConditions = CalcOrbitConditions(timeWait,order,direction,orbitalElement,uExhaust,asteroidCloseApp);
    else
        %Set asteroid conditions to close approach data
        asteroidConditions = asteroidCloseApp;
    end
    
    %Element is 1 to signify earth
    orbitalElement = 1;
    %Calculates conditions of earth at time of return to ESOI. The time
    %passed to it is the length of the heliocentric trip minus the wait
    %time.
    earthConditions = CalcOrbitConditions(timeGuess+timeWait,order,direction,orbitalElement,uExhaust,earthCloseApp);
    
elseif direction == Direction.FORWARD
    %Set asteroid conditions to close approach data
    asteroidConditions = asteroidCloseApp;
    %Element is 1 to signify earth
    orbitalElement = 1;
    %Calculates conditions of earth at time of escape. The time passed to
    %it is the length of the heliocentric trip.
    earthConditions = CalcOrbitConditions(timeGuess,order,direction,orbitalElement,uExhaust,earthCloseApp);
end

%A fourier series begins with 1 coefficient for order 0 and then adds 2 for
%every increase in order after that.  The last coefficient represents the
%earths escape angle from earth (toward sun, away from sun ect.) which is
%why there are 2*order+2 independent variables. The values being
%initialized to zero are arbitrary but we have found that they give the
%best results.

%Call fminSearchLanding in order to begin optimization
%Concatenating the variables that fminsearch will vary into 1 array
%Note this is a horizontal array
if optimize(LegNumber,1)
    x0 = [cConverge, departBound(LegNumber,1), approachBound(LegNumber,1), coastFraction(LegNumber,1)];
else
    x0 = cConverge;
end

[phiCoeff,fval,~,~,y,t] = fminSearchLanding(x0,order,direction,uExhaust,finalMass,efficiency,mdot,numberOfEngines,inputPower,escapeTime,timeGuess,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isSolar,COAST,departBound,approachBound,numIntervals,coastFraction,optimset('MaxFunEvals',Constants.FMIN_MAXNUMEVAL,'MaxIter',MaxIterations,'TolX',Constants.FMIN_TOLX,'TolFun',Constants.FMIN_TOLFUN),LegNumber,optimize,interpolateStruct,throttleFlag);
fprintf('Function   Value: %g\n',fval)

%Output information for user.  Also used for troubleshooting
%purposes.
if direction == Direction.RETURN
    TempF_R  = (y(end,1)-asteroidConditions(1))^2;
    TempF_Th = (y(end,2)-asteroidConditions(2))^2;
    TempF_U  = (y(end,3)-asteroidConditions(3))^2;
    TempF_H  = (y(end,4)-asteroidConditions(4))^2;
    fprintf('Function   Value Calculated from Y: %g\n',TempF_R+TempF_Th+TempF_U+TempF_H)
elseif direction == Direction.FORWARD
    RadiusSunSphere = earthConditions(1) + Constants.ESOI*cos(phiCoeff(end));
    forwardEscapeVelocity = sqrt((2*Constants.G * Constants.ME)/Constants.ESOI);
    TempF_R  = (y(end,1)-RadiusSunSphere)^2;
    TempF_Th = (y(end,2)-(earthConditions(2) + asin(sin(pi-phiCoeff(end))*Constants.ESOI/RadiusSunSphere)))^2;
    TempF_U  = (y(end,3)-(earthConditions(3) + forwardEscapeVelocity * sin(phiCoeff(end-1))))^2;
    TempF_H  = (y(end,4)-(Constants.EARTHH + Constants.ESOI * forwardEscapeVelocity * cos(phiCoeff(end-1))))^2;
    fprintf('Function   Value Calculated from Y: %g\n',TempF_R+TempF_Th+TempF_U+TempF_H)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Acceleration, Mass, and Thrust Outputs%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(1,1) = y(1,1);
% Y(2,1) = nan;
Y(2,1) = y(end,1);

if direction == Direction.RETURN
    finalMass  = finalMass + mdot*numberOfEngines*escapeTime*Constants.SCONVERSION;
end

for i=1:1:2
    
    
    %%%%%Set min force of thrust and efficiency
    minThrust = interpolateStruct(1).minThrust;
    minEfficiency = interpolateStruct(1).minEfficiency;
    
    
    
    powerInConstMDOT = interpolateStruct(1).powerInConstMDOT;
    etaConstMDOT =  interpolateStruct(1).etaConstMDOT;
    uExConstMDOT =  interpolateStruct(1).uExConstMDOT;
    mDotConstMDOT =  interpolateStruct(1).mDotConstMDOT;
    
    powerInConstUEX = interpolateStruct(1).powerInConstUEX;
    etaConstUEX =  interpolateStruct(1).etaConstUEX;
    uExConstUEX =  interpolateStruct(1).uExConstUEX;
    mDotConstUEX =  interpolateStruct(1).mDotConstUEX;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Find values %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %initialize so loop will begin
    conditionsCheck = 1;
    timeThrough = 0;
    conditionsNotMet = 0;
    
    while conditionsCheck > 0 %this loop will run no more than twice
        
        %check to see if Y is NAN
        if isnan(Y(i,1)) == 1
            eta = NaN;
            uExhaust = NaN;
            mdot = NaN;
            break;
        end
        
        
        
        
        %assume conditions will be met
        conditionsCheck = 0;
        
        %keep force of thrust higher (so keep MDOT constant)
        if throttleFlag == 0
            
            %%%%%%%%%%%%%%%%%%%%%%%Begin Index Search%%%%%%%%%%%%%%%%%
            
            %If the Power at r is larger than the power at the
            %earth(the highest tested value) then we set it to
            %the power at the earth. This is due to all the
            %thrusters and the solarpanels being constructed on
            %earth, therefore they can't be any more efficient
            %than when they are on earth.
            
            %determine power in based on r
            if y(i,1) < 1
                indexMDOT = 1;
                columnIndexMDOT = 1;
                %set Pin_r to be the max imput power from the
                %PowerVector
                Pin_r = powerInConstMDOT(indexMDOT, columnIndexMDOT);
                
            else
                %set Pin_r to be the max imput value from the
                %powervector divided by r squared
                Pin_r = powerInConstMDOT(1,1)/y(i,1)^2;
                %This is finding the columnIndex for constant
                %mdot.
                minPin = min(powerInConstMDOT,[],1);
                maxPin = max(powerInConstMDOT, [],1);
                for x= 1:1:size(powerInConstMDOT,2)
                    
                    if Pin_r >= minPin(1,x)
                        columnIndexMDOT = x;
                        break;
                    end
                end
                
                
                indexMDOT = (NUMBER_OF_INDICIES+1) - ((Pin_r/maxPin(1,columnIndexMDOT))*NUMBER_OF_INDICIES);
                
                
            end
            
            rowIndexMDOT = round(indexMDOT);
            
            
            %We do not expect this if to be necessary, but it is a
            %sanity check.
            if rowIndexMDOT < 1
                rowIndexMDOT = 1;
            elseif rowIndexMDOT > NUMBER_OF_INDICIES
                rowIndexMDOT = NUMBER_OF_INDICIES;
            end
            
            %%%%%%%%%%%%%%%%%%%%%End Index Search%%%%%%%%%%%%%%%%%%%%%
            
            %Use the index obtained to get uExhaust, mDot, eta for that instant
            uExhaust = uExConstMDOT(rowIndexMDOT, columnIndexMDOT);
            mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
            eta = etaConstMDOT(rowIndexMDOT,columnIndexMDOT);
            
            %Calculate thrust
            thrust = mdot * uExhaust;
            
            %if thrust or efficiecy is less than the min, move to
            %the next set of data
            while thrust < minThrust || eta < minEfficiency
                
                
                %if no values will fit the minimum requirements,
                %get out of loop and set a flag
                if (thrust < minThrust || eta < minEfficiency) && columnIndexMDOT == size(powerInConstMDOT,2)
                    
                    if timeThrough == 0 %first time through the loop
                        conditionsCheck = 1; %this will signal to loop back to the beggining
                        throttleFlag = 1; %change throttleFlag to use Const UEX data
                        
                        
                    elseif timeThrough == 1 %after changing data sets to find points that fit minimum conditions
                        
                        %if no values can be found that fit the minimum
                        %eta and thrust values return throttleFlag to
                        %its orgininal position to compare
                        throttleFlag = 0;
                        conditionsNotMet = 1;
                        
                    end
                    timeThrough = 1; %set time through to greater than 0
                    break;
                else
                    columnIndexMDOT = columnIndexMDOT +1;
                    indexMDOT= (NUMBER_OF_INDICIES+1) - ((Pin_r/maxPin(1,columnIndexMDOT))*NUMBER_OF_INDICIES);
                    rowIndexMDOT = round(indexMDOT);
                    
                    %We do not expect this if to be necessary, but it is a
                    %sanity check.
                    if rowIndexMDOT < 1
                        rowIndexMDOT = 1;
                    elseif rowIndexMDOT > NUMBER_OF_INDICIES
                        rowIndexMDOT = NUMBER_OF_INDICIES;
                    end
                    
                    uExhaust = uExConstMDOT(rowIndexMDOT, columnIndexMDOT);
                    mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
                    eta = etaConstMDOT(rowIndexMDOT,columnIndexMDOT);
                    
                    %Calculate thrust
                    thrust = mdot * uExhaust;
                end
                
                
                
            end
            
            
            %keep efficiency higher (so keep UEX as constant as possible)
        elseif throttleFlag == 1
            
            %%%%%%%%%%%%%%%%%Begin Index Search%%%%%%%%%%%%%%%%%
            
            %If the Power is larger than the largest tested
            %power then set it to the largest tested power
            if y(i,1) < 1
                indexUEX = 1;
                columnIndexUEX = 1;
                Pin_r = powerInConstUEX(indexUEX, columnIndexUEX);
            else
                %set Pin_r to be the max imput value from the
                %powervector divided by r squared
                Pin_r = powerInConstUEX(1,1)/y(i,1)^2;
                %Search for the Column index
                minPin = min(powerInConstUEX,[],1);
                maxPin = max(powerInConstUEX, [],1);
                for x= 1:1:size(powerInConstUEX,2)
                    
                    if Pin_r >= minPin(1,x)
                        columnIndexUEX = x;
                        break;
                    end
                end
                
                
                indexUEX = (NUMBER_OF_INDICIES+1) - ((Pin_r/maxPin(1,columnIndexUEX))*NUMBER_OF_INDICIES);
            end
            
            
            rowIndexUEX = round(indexUEX);
            %We do not expect this if to be necessary, but it is a
            %sanity check.
            if rowIndexUEX < 1
                rowIndexUEX = 1;
            elseif rowIndexUEX > NUMBER_OF_INDICIES
                rowIndexUEX = NUMBER_OF_INDICIES;
            end
            %%%%%%%%%%%%%%%%%%%%%End Index Search%%%%%%%%%%%%%%%%%%%%%
            
            %Use the index obtained to get uExhaust, mDot, eta for that instant
            uExhaust = uExConstUEX(rowIndexUEX, columnIndexUEX);
            mdot = mDotConstUEX(rowIndexUEX, columnIndexUEX);
            eta = etaConstUEX(rowIndexUEX,columnIndexUEX);
            
            %Calculate thrust
            thrust = mdot * uExhaust;
            
            while thrust < minThrust || eta < minEfficiency
                
                %if no values will fit the minimum requirements,
                %get out of loop and set a flag
                
                if (thrust < minThrust || eta < minEfficiency) && columnIndexUEX == size(powerInConstUEX,2)
                    if timeThrough == 0 %first time through the loop
                        conditionsCheck = 1; %this will signal to loop back to the beggining
                        throttleFlag = 0; %change throttleFlag to use Const UEX data
                        
                        
                    elseif timeThrough == 1 %after changing data sets to find points that fit minimum conditions
                        
                        %if no values can be found that fit the minimum
                        %eta and thrust values return throttleFlag to
                        %its orgininal position to compare
                        throttleFlag = 1;
                        conditionsNotMet = 1;
                        
                    end
                    timeThrough = 1; %set time through to greater than 0
                    
                    break;
                else
                    columnIndexUEX = columnIndexUEX +1;
                    indexUEX = (NUMBER_OF_INDICIES+1) - ((Pin_r/maxPin(1,columnIndexUEX))*NUMBER_OF_INDICIES);
                    rowIndexUEX = round(indexUEX);
                    
                    %We do not expect this if to be necessary, but it is a
                    %sanity check.
                    if rowIndexUEX < 1
                        rowIndexUEX = 1;
                    elseif rowIndexUEX > NUMBER_OF_INDICIES
                        rowIndexUEX = NUMBER_OF_INDICIES;
                    end
                    
                    uExhaust = uExConstUEX(rowIndexUEX, columnIndexUEX);
                    mdot = mDotConstUEX(rowIndexUEX, columnIndexUEX);
                    eta = etaConstUEX(rowIndexUEX,columnIndexUEX);
                    %Calculate thrust
                    thrust = mdot * uExhaust;
                    
                end
                
                
            end
        end
    end
    
    
    
    
    
    if conditionsNotMet == 1 %NO VALUES FOUND THAT MEET MINIMUM REQUIRMENTS
        %compare etas and thrust values for constant UEX and
        %MDOT, and choose values to use based on throttleFlag
        
        thrustUEX = uExConstUEX(rowIndexUEX, columnIndexUEX) * mDotConstUEX(rowIndexUEX, columnIndexUEX);
        etaUEX = etaConstUEX(rowIndexUEX, columnIndexUEX);
        thrustMDOT = uExConstMDOT(rowIndexMDOT, columnIndexMDOT) * mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
        etaMDOT = etaConstMDOT(rowIndexMDOT, columnIndexMDOT);
        
        if throttleFlag == 0
            if thrustUEX > thrustMDOT
                mdot = mDotConstUEX(rowIndexUEX, columnIndexUEX);
                uExhaust = uExConstUEX(rowIndexUEX, columnIndexUEX);
                eta = etaConstUEX(rowIndexUEX, columnIndexUEX);
            else
                mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
                uExhaust = uExConstantMDOT(rowIndexMDOT, columnIndexMDOT);
                eta = etaConstMDOT(rowIndexMDOT, columnIndexMDOT);
                
            end
        elseif throttleFlag ==1
            if etaUEX > etaMDOT
                mdot = mDotConstUEX(rowIndexUEX, columnIndexUEX);
                uExhaust = uExConstUEX(rowIndexUEX, columnIndexUEX);
                eta = etaConstUEX(rowIndexUEX, columnIndexUEX);
                
            else
                mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
                uExhaust = uExConstMDOT(rowIndexMDOT, columnIndexMDOT);
                eta = etaConstMDOT(rowIndexMDOT, columnIndexMDOT);
                
            end
        end
    end
    
    
    
    eff_r = eta;
    TempuExhaust = uExhaust;
    Tempmdot = mdot;
    
    % mdot calculation for constant uExhaust. Divided by
    % numberOfEngines in order to keep mdot as the mdot per
    % thruster
    
    %If the value for mdot is greater than is possible set it
    %at the maximum value
    
    if i==1
        fprintf('Final Acceleration = %g\n', (Tempmdot*TempuExhaust*numberOfEngines)/(finalMass - y(i,5)))%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Final num = %g\n', (Tempmdot*TempuExhaust*numberOfEngines))
        fprintf('Final denom = %g\n', (finalMass - y(i,5)))
    elseif i==2
        fprintf('Initial Acceleration = %g\n', (Tempmdot*TempuExhaust*numberOfEngines)/(finalMass - y(end,5)))
        fprintf('Initial num = %g\n', (Tempmdot*TempuExhaust*numberOfEngines))
        fprintf('Initial denom = %g\n', (finalMass - y(end,5)))
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%End Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sets 'interesting' based on whether the current value of fval returned by
%fminsearch is below the interesting threshold
if fval<Constants.FMIN_INTERESTINGTHRESHOLD
    isInteresting = 1;
else
    isInteresting = 0;
end

%If fminsearch evaluates an F which is less than G_fmin, the time is said
%to converge. This means that the the spacecraft would reach its
%destination. However, MinimizeTripTime will continue to run this function
%until it finds the lowest converging time guess

%Only use if fval<1 if trying to verfiy program works correctly, this will
%cause the craft to not land correctly if changed
%if fval<10
if FS== 1 %Full
    if fval<Constants.FMIN_FVALMAX
        %Update cConverge to the phi coefficients of the most recent
        %convergence.
        cConverge = phiCoeff(1,1:2*order+3);
        if COAST(LegNumber,1) && optimize(LegNumber,1)
            departBound(LegNumber,1) = phiCoeff(1,2*order+4);
            approachBound(LegNumber,1) = phiCoeff(1,2*order+5);
            coastFraction(LegNumber,1) = phiCoeff(1,2*order+6);
        end
        %Update other variables as a result of convergence
        converge = true;
        lastConverge = timeGuess;
        lastEscape = escapeTime;
        lastfval = fval;
    else
        converge = false;
    end
elseif FS==0
    if fval<Constants.FMIN_FVALMAXFAST
        %Update cConverge to the phi coefficients of the most recent
        %convergence.
        cConverge = phiCoeff(1,1:2*order+3);
        if COAST(LegNumber,1) && optimize(LegNumber,1)
            departBound(LegNumber,1) = phiCoeff(1,2*order+4);
            approachBound(LegNumber,1) = phiCoeff(1,2*order+5);
            coastFraction(LegNumber,1) = phiCoeff(1,2*order+6);
        end
        %Update other variables as a result of convergence
        converge = true;
        lastConverge = timeGuess;
        lastEscape = escapeTime;
        lastfval = fval;
    else
        converge = false;
    end
end

if converge==0
    fprintf('DID NOT CONVERGE\n\n')
elseif converge==1
    fprintf('CONVERGED\n\n')
end
end
