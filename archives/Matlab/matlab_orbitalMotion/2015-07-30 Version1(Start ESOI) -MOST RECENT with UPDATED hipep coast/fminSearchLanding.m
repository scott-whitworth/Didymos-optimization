%FMINSEARCHLANDING Calls fminsearch with landingObjectiveFunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function acts as a wrapper for landingObjectiveFunction and calls
%fminsearch with landingObjectiveFunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by fminsearch in timeObjectiveFunction
function [phiCoeff,fVal,exitFlag,output,MinYdata,MinTdata] = fminSearchLanding(x0,order,direction,uExhaust,finalMass,efficiency,mdot,numberOfEngines,inputPower,escapeTime,timeGuess,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,isSolar,COAST,passeddepartBound,passedapproachBound,numIntervals,passedcoastFraction,options,LegNumber,optimize,interpolateStruct,throttleFlag)

%These values are initialized for storage of the minimized values for the
%current time guess. They are initialized so that if it hits the max
%function evals it doesn't break the outputs in timeObjectiveFunction
lowestF = (rand(1)+.5)*1e20;
MinYdata = [nan,nan,nan,nan,nan];
MinTdata = nan;

%fminsearch is an unconstrained optamization algorythm in MATLAB.
[phiCoeff,fVal,exitFlag,output] = fminsearch(@landingObjectiveFunction,x0,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is passed to fminsearch and acts as the objective function
%and as such has to conform to a strict format where it takes a vector of
%independent variables and returns a funciton value. The purpose of the
%function is to evaluate the coefficients it is passed such that the global
%minimum of the function occurs when the ship lands on the asteroid. To do
%that it integrates the ships motion and then combines all of the ship and
%asteroid final positions in one number that represents error.

%Nan's are added to time and position data to serve as a 'divider' between
%coasting and thrusting segments.  This is used when T and Y are passed to
%Integrate for plot to use the segments correctly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by fminsearch in fminSearchLanding
    function F = landingObjectiveFunction(x0)
        
        integrationFlag = true;
        cConverge = x0(1,1:2*order+3);
        
        if COAST(LegNumber,1)
            if optimize(LegNumber,1)
                departBound = x0(1,2*order+4);
                approachBound = x0(1,2*order+5);
                coastFraction = x0(1,2*order+6);
                if departBound <= .05 || departBound >= 0.4 || approachBound <= .05 || approachBound >= .2 || coastFraction <= Constants.EPSILON || coastFraction >= 1-Constants.EPSILON
                    F_R  = nan;
                    F_Th = nan;
                    F_U  = nan;
                    F_H  = nan;
                    F = (rand(1)+0.5)*1e20;
                    integrationFlag = false;
                end
            else
                departBound = passeddepartBound(LegNumber,1);
                approachBound = passedapproachBound(LegNumber,1);
                coastFraction = passedcoastFraction(LegNumber,1);
            end
        else
            departBound = nan;
            approachBound = nan;
            coastFraction = nan;
        end
        
        %Function integrates the ships motion at a given time length for a
        %given set of phi coefficients (cConverge). It returns vectors for time and
        %position variables but only the position variables interest us.
        
        % For comments, see MinimizeTripTime
        % Note: integration bounds and code logic must be identical to that
        % in MinimizeTripTime and vice versa.
        if integrationFlag
            %Integrate non-coasting trips
            if ~COAST(LegNumber,1)
                initialPositionFlag = true;
                object = Object.SHIP_SUN_SPHERE;
                if direction == Direction.RETURN;
                    %Update ShipMass after escaping earth on the return
                    %trip.
                    ShipMass = finalMass + (mdot * numberOfEngines * escapeTime * Constants.SCONVERSION);
                    %Integrate the heliocentric portion of the return
                    [Tdata,Ydata] = IntegrateShipMotion(cConverge,[timeGuess,0],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,earthConditions,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                elseif direction == Direction.FORWARD
                    %Update shipmass to be the initial mass of the forward
                    %integration.  Since the integration is done in
                    %reverse, the final mass is used as the initial
                    %condition.
                    ShipMass = finalMass;
                    %Integrate the heliocentric portion of the forward
                    [Tdata,Ydata] = IntegrateShipMotion(cConverge,[timeGuess,0],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,asteroidConditions,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                end
                %Integrate Coasting Trips
            elseif COAST(LegNumber,1)
                %Determine time segments based on input/optimized
                %fractions.
                departThrustTime = departBound * timeGuess;
                approachThrustTime = approachBound * timeGuess;
                coastingSegment = timeGuess - departThrustTime - approachThrustTime;
                interval = coastingSegment/numIntervals(LegNumber,1);
                %Integrate return trip
                if direction == Direction.RETURN
                    %Set initial conditions
                    firstBound = timeGuess - approachThrustTime;
                    coastBound = firstBound - coastFraction * interval;
                    object = Object.SHIP_SUN_SPHERE;
                    initialPositionFlag = true;
                    ShipMass = finalMass + (mdot * numberOfEngines * escapeTime * Constants.SCONVERSION);
                    %Integrate approach
                    [firstT,firstY] = IntegrateShipMotion(cConverge,[timeGuess, firstBound],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,earthConditions,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                    Tdata = firstT;
                    %The final position of the spacecraft after this
                    %integration will be the initial conditions for the
                    %next segment
                    ShipPosition = firstY(end,:);
                    massExpended = firstY(end,5);
                    initialPositionFlag = false;
                    Ydata = firstY;
                    
                    %Loop through coasting segments, updating
                    %initial conditions and motion information at each
                    %iteration.
                    for timeSegment = 1:1:numIntervals(LegNumber,1)
                        %SecondY denotes a coasting segment, and to
                        %initiate coasting the object is changed to a
                        %NO_THRUST object
                        object = Object.NO_THRUST;
                        [secondT,secondY] = IntegrateShipMotion(cConverge,[firstBound, coastBound],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess, interpolateStruct,throttleFlag);
                        ShipPosition = horzcat(secondY(end,1:4),massExpended);
                        %ThirdY denotes a thrusting segment, and to thrust
                        %the object needs to be SHIP_SUN_SPHERE
                        object = Object.SHIP_SUN_SPHERE;
                        [thirdT,thirdY] = IntegrateShipMotion(cConverge,[coastBound, firstBound - interval],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                        %Time boundary variables are changed with each new
                        %coasting segment
                        firstBound = firstBound - interval;
                        coastBound = firstBound - coastFraction * interval;
                        ShipPosition = thirdY(end,:);
                        massExpended = thirdY(end,5);
                        %concatinate each segment of data into one variable
                        Ydata = vertcat(Ydata,[nan,nan,nan,nan,nan],secondY,[nan,nan,nan,nan,nan],thirdY);
                        Tdata = vertcat(Tdata,nan,secondT,nan,thirdT);
                        clear secondY;
                        clear thirdY;
                        clear secondT;
                        clear thirdT;
                    end
                    %Integrate departure from asteroid
                    [fourthT,fourthY] = IntegrateShipMotion(cConverge,[firstBound,0],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                    Ydata = vertcat(Ydata,[nan,nan,nan,nan,nan],fourthY);
                    Tdata = vertcat(Tdata,nan,fourthT);
                    
                    %See return for comments
                elseif direction == Direction.FORWARD
                    firstBound = timeGuess - approachThrustTime;
                    coastBound = firstBound - coastFraction * interval;
                    object = Object.SHIP_SUN_SPHERE;
                    initialPositionFlag = true;
                    ShipMass = finalMass;
                    [firstT,firstY] = IntegrateShipMotion(cConverge,[timeGuess, firstBound],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,asteroidConditions,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                    Tdata = firstT;
                    ShipPosition = firstY(end,:);
                    massExpended = firstY(end,5);
                    initialPositionFlag = false;
                    Ydata = firstY;
                    for timeSegment = 1:1:numIntervals(LegNumber,1)
                        object = Object.NO_THRUST;
                        [secondT,secondY] = IntegrateShipMotion(cConverge,[firstBound, coastBound],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                        ShipPosition = horzcat(secondY(end,1:4),massExpended);
                        object = Object.SHIP_SUN_SPHERE;
                        [thirdT,thirdY] = IntegrateShipMotion(cConverge,[coastBound, firstBound - interval],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                        firstBound = firstBound - interval;
                        coastBound = firstBound - coastFraction * interval;
                        ShipPosition = thirdY(end,:);
                        massExpended = thirdY(end,5);
                        Ydata = vertcat(Ydata,[nan,nan,nan,nan,nan],secondY,[nan,nan,nan,nan,nan],thirdY);
                        Tdata = vertcat(Tdata,nan,secondT,nan,thirdT);
                        clear secondY;
                        clear thirdY;
                        clear secondT;
                        clear thirdT;
                    end
                    [fourthT,fourthY] = IntegrateShipMotion(cConverge,[firstBound, 0],order,direction,uExhaust,ShipMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,ShipPosition,isSolar,object,initialPositionFlag,timeGuess,interpolateStruct,throttleFlag);
                    Ydata = vertcat(Ydata,[nan,nan,nan,nan,nan],fourthY);
                    Tdata = vertcat(Tdata,nan,fourthT);
                end
            end
            
            %This is the objective function. It is minimized when the final values of
            %radial position (r), angular postion (theta), radial velocity (u),
            %and specific angular momentum (h) are equal to those of the asteroid.
            %The accepted value of F is determined by G_fmin
            
            
            
            
            if direction == Direction.RETURN
                Fflag = true;
                %If the craft went within .5 AU (Ship burns up) it is
                %an unacceptable trajectory
                
                if min(Ydata(:,1)) < Constants.R_MIN
                    F_R  = nan;
                    F_Th = nan;
                    F_U  = nan;
                    F_H  = nan;
                    F = (rand(1)+0.5)*1e20;
                    Fflag = false;
                else
                    %Compare the position of the ship and the asteroid
                    F_R  = (Ydata(end,1)-asteroidConditions(1))^2;
                    F_Th = (Ydata(end,2)-asteroidConditions(2))^2;
                    F_U  = (Ydata(end,3)-asteroidConditions(3))^2;
                    F_H  = (Ydata(end,4)-asteroidConditions(4))^2;
                end
                
            else
                %For the Forward trip we are trying to hit ESOI instead of
                %the asteroid. Uses ESOI conditions based off of optimized
                %values for the escape angle (cConverge(end-1)) and the
                %angle of rotation of Earth's reference frame
                %(cConverge(end)).
                Fflag = true;
                forwardEscapeVelocity = sqrt((2*Constants.G * Constants.ME)/Constants.ESOI);
                if min(Ydata(:,1)) < Constants.R_MIN
                    F_R  = nan;
                    F_Th = nan;
                    F_U  = nan;
                    F_H  = nan;
                    F = (rand(1)+0.5)*1e20;
                    Fflag = false;
                else
                    %RadiusSunSphere is seperated as its own variable
                    %since it is used in several parameters
                    RadiusSunSphere = earthConditions(1) + Constants.ESOI*cos(cConverge(end));
                    F_R  = (Ydata(end,1)-RadiusSunSphere)^2;
                    F_Th = (Ydata(end,2)-(earthConditions(2) + asin(sin(pi-cConverge(end))*Constants.ESOI/RadiusSunSphere)))^2;
                    F_U  = (Ydata(end,3)-(earthConditions(3) + forwardEscapeVelocity * sin(cConverge(end-1))))^2;
                    F_H  = (Ydata(end,4)-(Constants.EARTHH + Constants.ESOI * forwardEscapeVelocity * cos(cConverge(end-1))))^2;
                end
                
            end
            if Fflag
                F = F_R + F_Th + F_U + F_H;
                %If the F value is the lowest F value so far store the
                %corresponding Ydata and Tdata vectors for plotting
                [lowestF,MinYdata,MinTdata] = StoreLowestFvalData(F,Ydata,Tdata,lowestF,MinYdata,MinTdata);
            end
        end
    end
end
