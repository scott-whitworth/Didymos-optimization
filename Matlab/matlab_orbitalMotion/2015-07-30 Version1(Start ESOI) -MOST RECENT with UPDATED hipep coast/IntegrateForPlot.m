%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the acceleration and mdot for the best Y and T data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FuelMass,lowestYdata,lowestTdata,Adata,departData,coastData,thrustData,approachData,MdotData,TmdotData] = IntegrateForPlot(COAST,direction,LegNumber,efficiency,mdot,numberOfEngines,inputPower,lowestEscape,lowestfinalMass,numIntervals,interpolateStruct,throttleFlag,lowestYdata,lowestTdata)

%Initialize coasting data
departData = [];
coastData = [];
thrustData = [];
approachData = [];

%If no coasting is desired for LegNumber, integrate from initial conditions
%for the entirety of lowestConverge.
if ~COAST(LegNumber,1)
    %For the non-coasting case the ship is considered to be
    %Object.SHIP_SUN_SPHERE for all acceleration calculations
    object = Object.SHIP_SUN_SPHERE;
    if direction == Direction.RETURN
        FuelMass = (mdot * numberOfEngines * (lowestEscape * Constants.SCONVERSION));
        %Calculate the mass of the ship at the end of the heliocentric
        %phase.  It is used as the final mass for accel calculations
        ShipMass = lowestfinalMass + FuelMass;
        %The fuel mass at the beginning of the trip is the sum from escape
        %and the helio portion of the trip
        FuelMass = FuelMass + abs(lowestYdata(end,5));
    elseif direction == Direction.FORWARD
        %The mass of the ship at the end of the helio portion of the trip
        %is just the final mass at the asteroid
        ShipMass = lowestfinalMass;
        %The fuel mass at the beginning of the trip is just the fuel mass
        %from the helio portion of the trip
        FuelMass = abs(lowestYdata(end,5));
    end
    
    %Acceleration cannot be passed back from ODE45 and fminsearch so it
    %must be calculated independently using the Y and T data passed
    %back for the lowest time
    Adata = GetAccel(lowestYdata,direction,object,ShipMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag);
    %Mdot must also be calculated independently. It is calculated from
    %the slope of y5
    [MdotData, TmdotData] = GetMdot(lowestYdata(:,5),lowestTdata,numberOfEngines,mdot);
    
    
    % If coasting is desired for LegNumber, coasting and non-coasting segments
    % must be integrated separately.
elseif COAST(LegNumber,1)
    if direction == Direction.RETURN
        FuelMass = (mdot * numberOfEngines * lowestEscape * Constants.SCONVERSION);
        ShipMass = lowestfinalMass + FuelMass;
        FuelMass = FuelMass + abs(lowestYdata(end,5));
    elseif direction == Direction.FORWARD
        ShipMass = lowestfinalMass;
        FuelMass = abs(lowestYdata(end,5));
    end
    
    %In order to calculate acceleration for coasting and seperate the
    %coasting and thrusting data the lowest Y vector needs to be seperated
    %back into it's components. The variable i is used to keep track of the
    %current index of lowestYdata while seperating the information.
    
    i=1;
    
    %NaNs are used in fminsearch are used to seperate each of the different
    %segments (firstY,secondY,thirdY and fourthY). To seperate the data it
    %loops until it finds a NaN saving the data to the appropriate vectors
    %along the way.
    while ~isnan(lowestYdata(i,1))
        approachData(i,:) = lowestYdata(i,:);
        i=i+1;
    end
    
    %For thrusting segments the object variable needs to be set to
    %Object.SHIP_SUN_SPHERE for proper acceleration calculations
    object = Object.SHIP_SUN_SPHERE;
    Adata = GetAccel(approachData,direction,object,ShipMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag);
    
    %The index i is increased by one to skip over the segment seperating
    %NaN value
    i=i+1;
    
    %It loops through all of the coasting and thrusting segments for the
    %current leg of the trip
    for timeSegment = 1:1:numIntervals(LegNumber,1)
        %The variable j is used to keep track of the current index for the
        %vectors containing information about individual segments and are
        %initialized to 1 before each new segment
        j=1;
        
        while ~isnan(lowestYdata(i,1))
            secondY(j,:) = lowestYdata(i,:);
            i=i+1;
            j=j+1;
        end
        
        %For coasting segments the object variable needs to be set to
        %Object.No_THRUST so get accel assigns a value of zero for
        %acceleration due to thrust
        object = Object.NO_THRUST;
        secondA = GetAccel(secondY,direction,object,ShipMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag);

        i=i+1;
        j=1;
        
        while ~isnan(lowestYdata(i,1))
            thirdY(j,:) = lowestYdata(i,:);
            i=i+1;
            j=j+1;
        end
        
        object = Object.SHIP_SUN_SPHERE;
        thirdA = GetAccel(thirdY,direction,object,ShipMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag);

        i=i+1;
        
        %This concatinates the acceleration vector with NaNs inbetween each
        %segment so that it's index matches that of Y and T
        Adata = vertcat(Adata,nan,secondA,nan,thirdA);
        
        %coastData and thrustData are seperated into seperate variables for
        %plotting the trip trajectory. Thrust data has a vector of NaNs
        %inbetween its values so that when the line is plotted it truncates
        %where coasting segments should be.
        coastData = vertcat(coastData,secondY);
        thrustData = vertcat(thrustData,[nan,nan,nan,nan,nan],thirdY);
        
        %The looping variables are cleared to prevent indexing errors 
        clear secondY;
        clear thirdY;
        clear thirdA;
        clear secondA;
    end
    
    j=1;
    
    sizeYdata = size(lowestYdata);
    %Since there is no NaN at the end of the lowestYdata vector the final
    %section is looped until i reaches the size of the vector
    while i <= sizeYdata(1,1)
        departData(j,:) = lowestYdata(i,:);
        i=i+1;
        j=j+1;
    end
    
    fourthA = GetAccel(departData,direction,object,ShipMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag);
    Adata = vertcat(Adata,nan,fourthA);
    
    %The data for mdot is calculated all at once
    [MdotData, TmdotData] = GetMdot(lowestYdata(:,5),lowestTdata,numberOfEngines,mdot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%