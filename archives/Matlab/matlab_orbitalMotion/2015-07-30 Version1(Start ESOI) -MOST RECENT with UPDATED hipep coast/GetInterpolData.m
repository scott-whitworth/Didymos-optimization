%This function reads in the Power versus uExhaust data from the
%excel files to be interpolated in solveOde45.  mdot can be
%calculated in solveode45 from the uExhaust data. (or vice versa, but we
%chose to do it this way)

function [interpolateStruct] = GetInterpolData(filename,thrusterNum,chooseExtrap,extrapPercent)


filename = strcat('InterpolData/',filename,'.txt');

%Read in interpolation data
InterpolData = fopen(filename);
temp = textscan(InterpolData,'%s','delimiter','\n');
%temp2 = textscan(InterpolData, '\', 'delimiter', '\n')
fclose(InterpolData);
RawInput = temp{1};


%With mDot constant
constantMDotVector = zeros(size(RawInput, 1)-1, 1);

% %Populate powerIn, mDot, Uex, and eta values for constant mDot

%potentially figure out how to declare RawInput2 size before loop occurs


for i=2:1:size(RawInput,1)
    
    RawInputSplit = strsplit(RawInput{i}, '\');
    if length(RawInputSplit) >1
        
        RawInput2(i,1) = RawInputSplit(1,2);
        
    end
    [temp,remainingTerm] = strtok(RawInput{i});
    constantMDotVector(i-1,1) = str2double(temp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantMDotVector(i-1,2) = str2double(utemp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantMDotVector(i-1,3) = str2double(utemp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantMDotVector(i-1,4) = str2double(utemp);
end




% With uExhaust constant

constantUexVector = zeros(size(RawInput2, 1)-1, 1);

%Populate powerIn, mDot, Uex, and eta values for constant uExaust
for i=2:1:size(RawInput2,1)
    [temp,remainingTerm] = strtok(RawInput2{i});
    constantUexVector(i-1,1) = str2double(temp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantUexVector(i-1,2) = str2double(utemp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantUexVector(i-1,3) = str2double(utemp);
    [utemp,remainingTerm] = strtok(remainingTerm);
    constantUexVector(i-1,4) = str2double(utemp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interplate constant Mdot vector


i = 1;
x = 1;
while constantMDotVector(i,1) ~= 0
    
    n = 1;
    while constantMDotVector(i,2) == constantMDotVector(i+1,2)
        
        powerIn(x,n) = constantMDotVector(i,1);
        mdot(x,n) = constantMDotVector(i,2);
        uExhaust(x,n) =  constantMDotVector(i,3);
        eta(x,n) = constantMDotVector(i,4);
        
        n = n + 1;
        i = i + 1;
    end
    
    powerIn(x,n) = constantMDotVector(i,1);
    mdot(x,n) = constantMDotVector(i,2);
    uExhaust(x,n) = constantMDotVector(i,3);
    eta(x,n) = constantMDotVector(i,4);
    
    x = x+1;
    i = i+1;
    
    if constantMDotVector(i,1) == 0
        break;
    end
end



%change all zeros to NaN

for i= 1:1:size(powerIn,2)
    for j = 1:1:size(powerIn,1)
        
        if powerIn(j,i) < 1
            powerIn(j,i) = NaN;
            
        end
    end
end


%set max and min for each set of constant mDot
minPower = min(powerIn,[],2);
maxPower = max(powerIn,[],2);



% if exprapolating set the maxPower to higher values
if chooseExtrap == 1
    maxPower = maxPower + (maxPower* extrapPercent) ;
end


%set the last minimum value to zero in order to interpolate
k = size(minPower,1);
minPower(k,1) = 0;

%Step for loop to calculate interpolation
stepSizePower = (maxPower-minPower)/Constants.NUMBER_OF_INDICIES;



%change all NaNs back to zero so Interp1 can function

for i= 1:1:size(powerIn,2)
    for j = 1:1:size(powerIn,1)
        a = isnan(powerIn);
        if a(j,i) == 1
            powerIn(j,i) = 0;
            
        end
    end
end




%set size of intererpolated power vector, eta, and mDot and Uex

PowerInterpMDOT = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn,1));
etaInterpMDOT = zeros(Constants.NUMBER_OF_INDICIES +1,size(powerIn,1));
etaInterp1 = zeros(Constants.NUMBER_OF_INDICIES ,1);
uExhaustInterpMDOT = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn,1));


%Set mDotInterp, each column has a different constant mDot
mDotInterpMDOT = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn,1));



mdotSize = size(mdot);
mdotCol = mdot(:,1)';
for i=1:size(PowerInterpMDOT,1)+1
    for j = 1:mdotSize(1,1)
        mDotInterpMDOT(i,j) = mdotCol(1,j);
    end
end


for j = 1: 1: size(powerIn,1)
    
    x = 1;
    for i = maxPower(j,1): -stepSizePower(j,1):minPower(j,1)
        PowerInterpMDOT(x,j) = i;
        
        x = x+1;
        
    end
    PowerValues = powerIn(j,:);
    etaValues = eta(j,:);
    PowerValues = PowerValues';
    etaValues = etaValues';
    
    a = length(PowerValues);
    for i=a:-1:1
        if PowerValues(i) <= 0
            PowerValues(i) =[];
            etaValues(i) = [];
        end
    end
    
    % put Zero as the minimum value before interpolating
    if j == size(powerIn,1)
        PowerValues(length(PowerValues) + 1,1) = 0;
        etaValues(length(etaValues) + 1,1) = 0;
        
    end
    %extrapolate eta values, term in '' is the type of MatLab extrapolation 
      if chooseExtrap == 1 
          etaInterp1 = interp1(PowerValues,etaValues,PowerInterpMDOT(:,j), 'pchip');
          
      else 
      etaInterp1 = interp1(PowerValues,etaValues,PowerInterpMDOT(:,j));
    
      end
      
    etaInterpMDOT(:,j) = etaInterp1;
    
end


% etaInterpMDOT



for i=1:1:size(PowerInterpMDOT,1)
    for j =1:1:size(PowerInterpMDOT,2)
        
        uExhaustInterpMDOT(i,j) = sqrt((2 * etaInterpMDOT(i,j) * PowerInterpMDOT(i,j))/ mDotInterpMDOT(i,j));
    end
end




% uExhaustInterpMDOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%% Interpolate for Constant uExhaust


i = 1;
x = 1;
while constantUexVector(i,1) ~= 0
    
    n = 1;
    while constantUexVector(i,3) == constantUexVector(i+1,3)
        
        powerIn2(x,n) = constantUexVector(i,1);
        mdot2(x,n) = constantUexVector(i,2);
        uExhaust2(x,n) =  constantUexVector(i,3);
        eta2(x,n) = constantUexVector(i,4);
        
        n = n + 1;
        i = i + 1;
    end
    
    powerIn2(x,n) = constantUexVector(i,1);
    mdot2(x,n) = constantUexVector(i,2);
    uExhaust2(x,n) = constantUexVector(i,3);
    eta2(x,n) = constantUexVector(i,4);
    
    x = x+1;
    i = i+1;
    
    if constantUexVector(i,1) == 0
        break;
    end
end




%change all zeros to NaN

for i= 1:1:size(powerIn2,2)
    for j = 1:1:size(powerIn2,1)
        
        if powerIn2(j,i) < 1
            powerIn2(j,i) = NaN;
            
        end
    end
end

%set max and min for each set of constant mDot
minPower = min(powerIn2,[],2);
maxPower = max(powerIn2,[],2);



% if exprapolating set the maxPower to higher values
if chooseExtrap == 1
    maxPower = maxPower + (maxPower* extrapPercent) ;
end



%set the last minimum value to zero in order to interpolate
k = size(minPower,1);
minPower(k,1) = 0;

%Step for loop to calculate interpolation
stepSizePower = (maxPower-minPower)/Constants.NUMBER_OF_INDICIES;



%change all NaNs back to zero so Interp1 can function

for i= 1:1:size(powerIn2,2)
    for j = 1:1:size(powerIn2,1)
        a = isnan(powerIn2);
        if a(j,i) == 1
            powerIn2(j,i) = 0;
            
        end
    end
end




%set size of intererpolated power vector, eta, and mDot and Uex

PowerInterpUX = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn2,1));
etaInterpUX = zeros(Constants.NUMBER_OF_INDICIES +1,size(powerIn2,1));
etaInterp1 = zeros(Constants.NUMBER_OF_INDICIES ,1);
mDotInterpUX = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn2,1));


%Set uExhaustInterp, each column has a different constant uExhaust

uExhaustInterpUX = zeros(Constants.NUMBER_OF_INDICIES,size(powerIn2,1));

uExSize = size(uExhaust2);
uExhaustCol = uExhaust2(:,1)';
for i=1:size(PowerInterpUX,1)+1
    for j = 1:uExSize(1,1)
        uExhaustInterpUX(i,j) = uExhaustCol(1,j);
    end
end


for j = 1: 1: size(powerIn2,1)
    
    x = 1;
    for i = maxPower(j,1): -stepSizePower(j,1):minPower(j,1)
        PowerInterpUX(x,j) = i;
        
        x = x+1;
        
    end
    PowerValues = powerIn2(j,:);
    etaValues = eta2(j,:);
    PowerValues = PowerValues';
    etaValues = etaValues';
    
    
    
    a = length(PowerValues);
    for i=a:-1:1
        if PowerValues(i) == 0
            PowerValues(i) =[];
            etaValues(i)= [];
        end
    end
    
    % put Zero as the minimum value before interpolating
    if j == size(powerIn2,1)
        PowerValues(length(PowerValues) + 1,1) = 0;
        etaValues(length(etaValues) + 1,1) = 0;
        
    end
    
   %extrapolate eta values, term in '' is the type of MatLab extrapolation 
      if chooseExtrap == 1 
      etaInterp1 = interp1(PowerValues,etaValues,PowerInterpUX(:,j), 'pchip');
          
      else 
      etaInterp1 = interp1(PowerValues,etaValues,PowerInterpUX(:,j));
    
      end
          
    etaInterpUX(:,j) = etaInterp1;
end





for i=1:1:size(PowerInterpUX,1)
    for j =1:1:size(PowerInterpUX,2)
        
        mDotInterpUX(i,j) = (2* etaInterpUX(i,j) * PowerInterpUX(i,j))/ (uExhaustInterpUX(i,j)^2);
        
    end
end


%plot(PowerInterpUX(:,1), etaInterpUX(:,1), 'r+',PowerInterpMDOT(:,1), etaInterpMDOT(:,1), 'bo');



%declare struct variable
interpolateStruct = struct();

interpolateStruct(1).powerInConstUEX = PowerInterpUX;
interpolateStruct(1).etaConstUEX = etaInterpUX;
interpolateStruct(1).uExConstUEX = uExhaustInterpUX;
interpolateStruct(1).mDotConstUEX = mDotInterpUX;

interpolateStruct(1).powerInConstMDOT = PowerInterpMDOT;
interpolateStruct(1).etaConstMDOT = etaInterpMDOT;
interpolateStruct(1).uExConstMDOT = uExhaustInterpMDOT;
interpolateStruct(1).mDotConstMDOT = mDotInterpMDOT;


interpolateStruct(1).thusterID = thrusterNum;
interpolateStruct(1).minThrust = xlsread('Input.xlsx','Operation Guidelines','b1:b1');
interpolateStruct(1).minEfficiency = xlsread('Input.xlsx','Operation Guidelines','b2:b2');





end


