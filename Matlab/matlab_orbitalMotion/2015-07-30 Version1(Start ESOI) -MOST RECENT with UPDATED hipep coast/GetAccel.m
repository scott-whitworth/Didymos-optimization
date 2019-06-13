function a = GetAccel(y,direction,object,finalMass,efficiency,mdot,numberOfEngines,inputPower,interpolateStruct,throttleFlag)
EffBoost = Constants.EFFBOOST;
numPoints = size(y);
NUMBER_OF_INDICIES = Constants.NUMBER_OF_INDICIES;

a = zeros(numPoints(1,1),1);

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





MaxMdot = mdot;
for i = 1:1:numPoints(1,1)
    if object == Object.NO_THRUST
        if direction == Direction.FORWARD
            a(i,1) = 0;
        else
            a(i,1) = 0;
        end
    elseif object == Object.SHIP_SUN_SPHERE
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Find values %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %initialize so loop will begin
        conditionsCheck = 1;
        timeThrough = 0;
        conditionsNotMet = 0;
        
        while conditionsCheck > 0 %this loop will run no more than twice
            
            %check to see if Y is NAN
            if isnan(y(i,1)) == 1
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
        
        %
        %         %If the value for mdot is greater than is possible set it
        %         %at the maximum value
        %         if mdot > MaxMdot
        %             mdot = MaxMdot;
        %         end
        
        % Generalized acceleration equation for any value of mdot
        % or uExhaust
        a(i,1) = (mdot*uExhaust*numberOfEngines)/(finalMass - y(i,5));
    end
end
end