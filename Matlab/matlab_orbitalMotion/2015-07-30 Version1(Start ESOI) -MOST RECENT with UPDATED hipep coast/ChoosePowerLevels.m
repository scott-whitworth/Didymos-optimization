%Allow user to choose which power levels to evaluate

%Called by MultipleRuns
function [powerStart,powerEnd] = ChoosePowerLevels()
exitFlag = false;
while ~exitFlag
    exitFlag = true;
    lowerpowerChoice = input('Lower engine multiplier: ');
    upperpowerChoice = input('Upper engine multiplier: ');
    if lowerpowerChoice <= 0 || upperpowerChoice <= 0 || mod(lowerpowerChoice,1) ~= 0 || mod(upperpowerChoice,1) ~= 0
        fprintf('\nFatal Inputs. Please enter integers greater than 0.\n')
        exitFlag = false;
    end
end
%Set start and end to same value in order iterate over one power level
powerStart = lowerpowerChoice;
powerEnd = upperpowerChoice;

end

