% Allow user to choose the Alpha P from a preset list.
% Currently set to prevent looping for different Alpha P values.

% Called by MultipleRuns
function [minAlphaP,maxAlphaP] = ChooseAlphaP()

exitFlag = false;
while ~exitFlag
    exitFlag = true;
    fprintf('\nSelect Alpha P value\n')
    fprintf('1  16\n')
    fprintf('2  32\n')
    fprintf('3  48\n')
    fprintf('4  16 & 32\n\n')
    
    choice = input(': ');
    if choice < 0 || choice > 4 || mod(choice,1) ~= 0
        fprintf('Fatal Input. Please enter an integer from the list.\n');
        exitFlag = false;
    end
end
switch choice
    case 1
        minAlphaP = 16;
        maxAlphaP = 16;
    case 2
        minAlphaP = 32;
        maxAlphaP = 32;
    case 3
        minAlphaP = 48;
        maxAlphaP = 48;
    case 4
        minAlphaP = 16;
        maxAlphaP = 32;
end