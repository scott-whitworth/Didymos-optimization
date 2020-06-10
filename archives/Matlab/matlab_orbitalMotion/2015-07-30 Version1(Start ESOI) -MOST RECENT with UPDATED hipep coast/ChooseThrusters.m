%Allow user to choose which thrusters to evaluate

%Called by MultipleRuns
function [numChoiceThrusters,thrusterStructsArr,nameofThruster] = ChooseThrusters(dataEntry, throttleFlag)

%Select from existing thrusters or input new thruster
if dataEntry == 0
    numChoiceThrusters = input('Input number of thrusters to be evaluated, input 0 for all thrusters, or -1 to add a thruster): ');
else
    numChoiceThrusters = xlsread('Input.xlsx','Base Inputs','b10:b10');
end
%Add a new thruster
if (numChoiceThrusters == -1)
    %User input thruster parameters
    thrusterName = input('Input thruster name: ','s');
    alphaT = input('Input AlphaT: ','s');
    numEnginesMultiplier = input('Input the multiplier for number of engines: ','s');
    inputPower = input('Input the input power for one engine: ','s');
    uex = input('Input Uex: ','s');
    efficiency = input('Input efficiency: ','s');
    
    %Open file for appending
   
  
    thrusterData = fopen('Thrusters.txt', 'a');
    
    
    %Print thruster data into file
    fprintf(thrusterData, '%s\t%s\t%s\t%s\t%s\t%s\r\n',thrusterName,alphaT,numEnginesMultiplier,inputPower,uex,efficiency);
    fclose(thrusterData);
end

%Open thruster data file
 %When using CONST MDOT and CONST UEX values, there are different max
    %values for each. Depending on the throttleFlag, (0 = constant MDOT and
    %1 = constant UEX) the code will call the corresponding text file. 
    if throttleFlag == 0 %constant MDOT
         
        thrusterData = fopen('ThrustersMDOT.txt');
    elseif throttleFlag == 1 %constant UEX
         thrusterData = fopen('ThrustersUEX.txt');
    end
    


%Reads the text file one line at a time, keeping the whole line as a
%string. 'temp' is a dummy cell array--it is not actually usable.
temp = textscan(thrusterData, '%s','delimiter','\n');
fclose(thrusterData);

%Reads values from 'temp' into a usable cell array 'RawInput' containing
%each line of the original data file as a string
RawInput = temp{1};

%Declare struct variable
thrusterStructs = struct();

%Convert text from file to thruster structs
for i = 1:1:size(RawInput,1)
    [thrusterName,remainingTerm] = strtok(RawInput{i});
    [alphaT,remainingTerm] = strtok(remainingTerm);
    alphaT = str2double(alphaT);
    [numEnginesMultiplier,remainingTerm] = strtok(remainingTerm);
    numEnginesMultiplier = str2double(numEnginesMultiplier);
    [inputPower,remainingTerm] = strtok(remainingTerm);
    inputPower = str2double(inputPower);
    [uex,remainingTerm] = strtok(remainingTerm);
    uex = str2double(uex);
    [efficiencyTemp,remainingTerm] = strtok(remainingTerm);
    efficiency(1) = str2double(efficiencyTemp);
    [efficiencyTemp,remainingTerm] = strtok(remainingTerm);
    efficiency(2) = str2double(efficiencyTemp);
    [efficiencyTemp,remainingTerm] = strtok(remainingTerm);
    efficiency(3) = str2double(efficiencyTemp);
    [efficiencyTemp,remainingTerm] = strtok(remainingTerm);
    efficiency(4) = str2double(efficiencyTemp);
    [efficiencyTemp,remainingTerm] = strtok(remainingTerm);
    efficiency(5) = str2double(efficiencyTemp);
    [heliotripfilename,remainingTerm] = strtok(remainingTerm);
    [earthtripfilename,remainingTerm] = strtok(remainingTerm);
    [helioplotfilename,remainingTerm] = strtok(remainingTerm);
    [earthplotfilename,remainingTerm] = strtok(remainingTerm);
    [interpoldatafilename,remainingTerm] = strtok(remainingTerm);
    thrusterStructs(i).ThrusterName = thrusterName;
    thrusterStructs(i).AlphaT = alphaT;
    thrusterStructs(i).NumEnginesMultiplier = numEnginesMultiplier;
    thrusterStructs(i).InputPower = inputPower;
    thrusterStructs(i).Uex = uex;
    thrusterStructs(i).Efficiency = efficiency;
    thrusterStructs(i).HelioTripFileName = heliotripfilename;
    thrusterStructs(i).EarthTripFileName = earthtripfilename;
    thrusterStructs(i).HelioPlotFileName = helioplotfilename;
    thrusterStructs(i).EarthPlotFileName = earthplotfilename;
    thrusterStructs(i).InterpolDataFileName = interpoldatafilename;
    thrusterStructs(i).ThrusterIDNumber = i;
end





%Holds the total number of thrusters
numThrusters = length(thrusterStructs);
%If chosen to evaluate certain thrusters
if (numChoiceThrusters > 0)
    thrusterChoice(numChoiceThrusters) = 0;
    if dataEntry == 0
        for i = 1:1:numThrusters
            fprintf('%i',i);
            fprintf([' for ' thrusterStructs(i).ThrusterName '\n']);
        end
    end
    
    %Step through number of thrusters to allow user to choose thrusters
    if dataEntry == 0
        for i = 1:1:numChoiceThrusters
            thrusterChoice(i) = input('Input thruster ID: ');
        end
    else
        thrusterChoice = xlsread('Input.xlsx','Base Inputs','b11:b21');
    end
    for i = 1:1:numChoiceThrusters
        nameofThruster{i} = thrusterStructs(thrusterChoice(i)).ThrusterName;
    end
    
    %Store selected thrusters in array to be returned
    for i = 1:1:numChoiceThrusters
        thrusterStructsArr(i) = thrusterStructs(thrusterChoice(i));
    end
    
    %If chosen to evaluate all thrusters
elseif (numChoiceThrusters == 0)
    
    %Set the choice of number of thrusters to the total number of thrusters
    numChoiceThrusters = numThrusters;
    thrusterStructsArr = thrusterStructs;
end


end

