% Call ode45 with ode45Fun to integrate the equations of motion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function acts as a wrapper for ode45Fun and calls ode45 with ode45Fun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by ode45 in EscapeEarth,EscapeEarthForward,CalcOrbitConditions, and IntegrateShipMotion
function [t,Y,te,xe,ie] = solveOde45Fun(time,initial,order,direction,object,phiCoeff,options,uExhaust,accelFinal,finalMass,efficiency,mdot,numberOfEngines,inputPower,isSolar,throttleFlag,minThrust,minEfficiency,powerInConstMDOT,etaConstMDOT,uExConstMDOT,mDotConstMDOT,powerInConstUEX,etaConstUEX,uExConstUEX,mDotConstUEX)

%MaxMdot = mdot;
%We commented MaxMdot out because we dont seem to be using
%it in the code. it is connected to line 62 and seems to be self contained.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ode45Dun is called vastly more often than any other function in
%the program, and thus it needs to run very efficiently. In order
%to accomplish this, calls to the classes Direction, Object, and
%Constants have been replaced with local definitions of these
%classes.
RETURN = Direction.RETURN;
EffBoost = Constants.EFFBOOST;
NO_THRUST = Object.NO_THRUST;
SHIP_EARTH_SPHERE = Object.SHIP_EARTH_SPHERE;
SHIP_SUN_SPHERE = Object.SHIP_SUN_SPHERE;
NUMBER_OF_INDICIES = Constants.NUMBER_OF_INDICIES;
MUS = Constants.MUS;
MUE = Constants.MUE;
MCONVERSION = Constants.MCONVERSION;
SCONVERSION = Constants.SCONVERSION;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('options', 'var'))
    options = [];
end

if nargout == 5
    [t,Y,te,xe,ie] = ode45(@ode45Fun, time, initial, options);
elseif nargout == 4
    [t,Y,te,xe] = ode45(@ode45Fun, time, initial, options);
elseif nargout == 3
    [t,Y,te] = ode45(@ode45Fun, time, initial, options);
elseif nargout == 2
    [t,Y] = ode45(@ode45Fun, time, initial, options);
elseif nargout == 1
    [t] = ode45(@ode45Fun, time, initial, options);
else
    ode45(@ode45Fun, time, initial, options);
end

%ODE45FUN Used by ode45 to integrate the equations of motion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is used by ode45 to integrate the equations of motion. It is
%called to calculate the derivatives of radial position (r), angular
%postion (theta), radial velocity (u), and specific angular momentum (h).
%It accesses G_object to know if it is working in sun or earth sphere of
%influence and to know if there is thrust and how to calculate that thrust.
%The indicator is set in functions outside this one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by solveOde45Fun
    function dy = ode45Fun(t,y)
        
        % mdot = MaxMdot;
        
        %Sets mu to correspond to the mass of the central body
        if object == NO_THRUST
            mu = MUS;
            accel = 0;
            phi = 0;
            mdot = 0;
            numberOfEngines = 0;
        elseif object == SHIP_EARTH_SPHERE
            mu = MUE;
            %Since the engine never shuts off, acceleration can be
            %calculated directly from current trip time. As fuel burns,
            %weight decreases and acceleration increases. Time is
            %maintained continuous throughout the program so it will always
            %represent current total trip time. For the return trip time is
            %running backward so weight will be increasing and acceleration
            %decreasing options for forward or return trip.
            if direction == RETURN
                %Must start with dryMass since integration begins at the
                %end of the actual trip. When EscapeEarth calls
                %solveOde45Fun, time moves in the positive direction for
                %both forward and return trips. Acceleration decreases in
                %computational time, but increases in real time.
                accel = accelFinal/(1+(accelFinal*(t*SCONVERSION)/uExhaust));
            end
            %Makes phi tangent to motion
            phi = atan(y(1)*y(3)/y(4));
            
        elseif object == SHIP_SUN_SPHERE
            mu = MUS;
            if isSolar == 1 %This is assuming that the solar power varies as 1/r^2.
                
                %The set of code below is duplicated in GetAccel and in
                %TimeObjectiveFun. If you are going to make a change to
                %either you must make a change to this set of code. This
                %cope of the code will have the most commenting.
                
                
                
                
                
                
                %             %%%%%Set min force of thrust and efficiency
                %             minThrust = interpolateStruct(1).minThrust;
                %             minEfficiency = interpolateStruct(1).minEfficiency;
                %
                %             %Struct created in GetInterpolData
                %             powerInConstMDOT = interpolateStruct(1).powerInConstMDOT;
                %             etaConstMDOT =  interpolateStruct(1).etaConstMDOT;
                %             uExConstMDOT =  interpolateStruct(1).uExConstMDOT;
                %             mDotConstMDOT =  interpolateStruct(1).mDotConstMDOT;
                %
                %             powerInConstUEX = interpolateStruct(1).powerInConstUEX;
                %             etaConstUEX =  interpolateStruct(1).etaConstUEX;
                %             uExConstUEX =  interpolateStruct(1).uExConstUEX;
                %             mDotConstUEX =  interpolateStruct(1).mDotConstUEX;
                %
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%% Find values %%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %initialize so loop will begin
                conditionsCheck = 1;
                timeThrough = 0;
                conditionsNotMet = 0;
                
                while conditionsCheck > 0 %this loop will run no more than twice
                    
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
                        if y(1) < 1
                            indexMDOT = 1;
                            columnIndexMDOT = 1;
                            %set Pin_r to be the max imput power from the
                            %PowerVector
                            Pin_r = powerInConstMDOT(indexMDOT, columnIndexMDOT);
                            
                        else
                            %set Pin_r to be the max imput value from the
                            %powervector divided by r squared
                            Pin_r = powerInConstMDOT(1,1)/y(1)^2;
                            %This is finding the columnIndex for constant
                            %mdot.
                            minPin = min(powerInConstMDOT,[],1);
                            maxPin = max(powerInConstMDOT, [],1);
                            for i= 1:1:size(powerInConstMDOT,2)
                                
                                if Pin_r >= minPin(1,i)
                                    columnIndexMDOT = i;
                                    break;
                                end
                            end
                            
                            %using the columnIndex and the maxPin(1,1) we
                            %find the rowIndex.
                            indexMDOT = (NUMBER_OF_INDICIES+1) - ((Pin_r/maxPin(1,columnIndexMDOT))*NUMBER_OF_INDICIES);
                            
                            
                        end
                        
                        rowIndexMDOT = round(indexMDOT);
                        %The row and column Indicies will point to
                        %equivalent numbers in the power matrix,
                        %mDot,uExhaust, and Eta matricies. All four numbers
                        %are one set.
                        
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
                                
                                %Using the row and column Index, the
                                %code finds the uExhaust,mDot and eta for
                                %later calculations.
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
                        if y(1) < 1
                            indexUEX = 1;
                            columnIndexUEX = 1;
                            Pin_r = powerInConstUEX(indexUEX, columnIndexUEX);
                        else
                            %set Pin_r to be the max imput value from the
                            %powervector divided by r squared
                            Pin_r = powerInConstUEX(1,1)/y(1)^2;
                            %Search for the Column index
                            minPin = min(powerInConstUEX,[],1);
                            maxPin = max(powerInConstUEX, [],1);
                            for i= 1:1:size(powerInConstUEX,2)
                                
                                if Pin_r >= minPin(1,i)
                                    columnIndexUEX = i;
                                    break;
                                end
                            end
                            
                            
                            %Using the column Index we find the row index
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
                            
                            if (thrust < minThrust || eta < minEfficiency)&& columnIndexUEX == size(powerInConstUEX,2)
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
                        else
                            mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
                            uExhaust = uExConstMDOT(rowIndexMDOT, columnIndexMDOT);
                        end
                    elseif throttleFlag ==1
                        if etaUEX > etaMDOT
                            mdot = mDotConstUEX(rowIndexUEX, columnIndexUEX);
                            uExhaust = uExConstUEX(rowIndexUEX, columnIndexUEX);
                        else
                            mdot = mDotConstMDOT(rowIndexMDOT, columnIndexMDOT);
                            uExhaust = uExConstMDOT(rowIndexMDOT, columnIndexMDOT);
                        end
                    end
                end
                
                % Generalized acceleration equation for any value of mdot
                % or uExhaust
                accel = (mdot*uExhaust*numberOfEngines)/(finalMass - y(5));
                
                
            else %This is assuming a constant power source
                %%%%%%%%%%%%%%%%%%IGNORE (Always Solar)%%%%%%%%%%%%%%%%%%
                %Must start with dryMass since integration begins at the
                %end of the actual trip. When fminSearchLanding calls
                %IntegrateShipMotion, negative direction for the return
                %trip, and in the positive direction for the forward trip.
                %Acceleration increases in real time.
                accel = accelFinal/(1-(accelFinal*(t*SCONVERSION)/uExhaust));
            end
            
            
            %Evaluates value of Fourier series
            %The first term is simply a0/2
            phi=phiCoeff(1)/2;
            for j=1:order
                %The following follows the pattern of the Fourier series
                %an*cos(nx)+bn*sin(nx) where the pattern in phi goes
                %firstvalue,a1,b1,a2,b2....an,bn,escape angle.
                phi=phi+phiCoeff(2*j)*cos(j*t/2)+phiCoeff(2*j+1)*sin(j*t/2);
                %%%%%%%%%%%%%%%%End IGNORE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        %Acceleration needs to be nondimensionalized
        accel = accel * SCONVERSION^2/MCONVERSION;
        %mdot needs to be nondimensionalized
        mdot = mdot * SCONVERSION;
        
        %Equations of orbital motion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initialize array
        dy = zeros(5,1);
        %Derivative of radial postion
        dy(1) = y(3);
        %Derivative of angular postion is directly related to r and h.
        dy(2) = y(4)/y(1)^2;
        %Derivative of radial velocity is composed of centrifugal force
        %gravity term + thrust term.
        dy(3) = y(4)^2/y(1)^3-(mu/y(1)^2)+ accel*sin(phi);
        %Derivative of specific angular momentum is only affected by thrust
        dy(4) = y(1)*accel*cos(phi);
        %Derivative of mass expended.
        dy(5) = mdot * numberOfEngines;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

end
