% Allow user to specify coasting parameters for both forward and return
% trips.

% Called by MultipleRuns
function [COAST,departBound,approachBound,numIntervals,coastFraction,optimize] = ChooseCoast(dataEntry)

exitFlag = false;
badexcelinputs = 0;
while ~exitFlag
    exitFlag = true;
    if dataEntry == 0
        fprintf('Would you like to coast?\n\n')
        fprintf('0 No\n')
        fprintf('1 Yes\n\n')
        choice = input(': ');
    else
        if badexcelinputs == 1
            error('Problem with coasting parameters in excel file.  Please verify that excel all numbers of in the correct range')
        end
        choice = xlsread('Input.xlsx','Operation Guidelines','b3:b3');
    end
    if choice < 0 || choice > 1 || mod(choice,1) ~= 0
        fprintf('Fatal Input. Please enter an integer from the list.\n');
        exitFlag = false;
        badexcelinputs = 1;
    end
end
% Initialize parameters
% First row holds return values because it is calculated first
% Second row holds forward values
% Designed so LegNumber can be passed to the vector from MultipleRuns
departBound = zeros(2,1);
approachBound = zeros(2,1);
numIntervals = zeros(2,1);
coastFraction = zeros(2,1);
COAST = zeros(2,1);
optimize = zeros(2,1);

switch choice
    % If no coasting is desired,set all coasting flags to false and all
    % coasting parameters to nan because they won't be used.
    case 0
        COAST(:,1) = false;
        departBound(:,1) = nan;
        approachBound(:,1) = nan;
        numIntervals(:,1) = nan;
        coastFraction(:,1) = nan;
        % If coasting is desired, then the user has the choice whether to coast
        % on the forward trip, return trip, or both.
    case 1
        badexcelinputs = 0;
        exitFlag = false;
        while ~exitFlag
            exitFlag = true;
            if dataEntry == 0
                fprintf('Would you like to coast during the forward trip?\n')
                fprintf('0 No\n')
                fprintf('1 Yes\n\n')
                forwardchoice = input(': ');
                if forwardchoice == 1
                    optimize(2,1) = input('Would you like to optimize parameters? No(0), Yes(1)\n');
                else
                    optimize(2,1) = 0;
                end
            else
                if badexcelinputs == 1
                    error('Problem with coasting parameters in excel file.  Please verify that excel all numbers of in the correct range')
                end
                forwardchoice = xlsread('Input.xlsx','Operation Guidelines','b4:b4');
                optimize(2,1) = xlsread('Input.xlsx','Operation Guidelines','b5:b5');
            end
            if forwardchoice < 0 || forwardchoice > 1 || mod(forwardchoice,1) ~= 0
                fprintf('Fatal Input. Please enter an integer from the list.\n');
                exitFlag = false;
                badexcelinputs = 1;
            end
            if optimize(2,1) < 0 || optimize(2,1) > 1 || mod(optimize(2,1),1) ~= 0
                fprintf('Fatal Input. Please enter an integer from the list.\n');
                exitFlag = false;
                badexcelinputs = 1;
            end
        end
        switch forwardchoice
            % If no coasting is desired for the forward trip,set forward
            % coasting flag to false and all forward coasting parameters
            % to nan because they won't be used.
            case 0
                COAST(2,1) = false;
                departBound(2,1) = nan;
                approachBound(2,1) = nan;
                numIntervals(2,1) = nan;
                coastFraction(2,1) = nan;
                % Allow user to choose forward coasting parameters
            case 1
                COAST(2,1) = true;
                exitFlag = false;
                badexcelinputs = 0;
                if ~optimize(2,1)
                    while ~exitFlag
                        exitFlag = true;
                        if dataEntry == 0
                            departBound(2,1) = input('\nEnter thrust time for departure (as a fraction of one-way trip time): ');
                            approachBound(2,1) = input('Enter thrust time for approach  (as a fraction of one-way trip time): ');
                            coastFraction(2,1) = input('Enter coasting fraction per interval forward: ');
                            numIntervals(2,1) = input('Enter number of coasting intervals forward: ');
                        else
                            if badexcelinputs == 1
                                error('Problem with coasting parameters in excel file.  Please verify that all numbers are in the correct range')
                            end
                            departBound(2,1) = xlsread('Input.xlsx','Operation Guidelines','b6:b6');
                            approachBound(2,1) = xlsread('Input.xlsx','Operation Guidelines','b7:b7');
                            coastFraction(2,1) = xlsread('Input.xlsx','Operation Guidelines','b9:b9');
                            numIntervals(2,1) = xlsread('Input.xlsx','Operation Guidelines','b8:b8');
                        end
                        if departBound(2,1) <= 0 || approachBound(2,1) <= 0 || departBound(2,1) >= 1 || approachBound(2,1) >= 1 || coastFraction(2,1) <= 0 || coastFraction(2,1) >= 1 || numIntervals(2,1) <= 0 || mod(numIntervals(2,1),1) ~= 0
                            fprintf('\nFatal Inputs. Please enter bounds between 0 and 1\n')
                            exitFlag = false;
                            badexcelinputs = 1;
                        end
                    end
                else
                    numIntervals(2,1) = xlsread('Input.xlsx','Operation Guidelines','b8:b8');
                    departBound(2,1) = nan;
                    approachBound(2,1) = nan;
                    coastFraction(2,1) = nan;
                end
        end
        badexcelinputs = 0;
        exitFlag = false;
        while ~exitFlag
            exitFlag = true;
            if dataEntry == 0
                fprintf('Would you like to coast during the return trip?\n')
                fprintf('0 No\n')
                fprintf('1 Yes\n\n')
                returnchoice = input(': ');
                if returnchoice == 1
                    optimize(1,1) = input('Would you like to optimize parameters? No(0), Yes(1)\n');
                else
                    optimize(1,1) = 0;
                end
            else
                if badexcelinputs == 1
                    error('Problem with coasting parameters in excel file.  Please verify that excel all numbers of in the correct range')
                end
                returnchoice = xlsread('Input.xlsx','Operation Guidelines','b10:b10');
                optimize(1,1) = xlsread('Input.xlsx','Operation Guidelines','b11:b11');
            end
            if returnchoice < 0 || returnchoice > 1 || mod(returnchoice,1) ~= 0
                fprintf('Fatal Input. Please enter an integer from the list.\n');
                exitFlag = false;
                badexcelinputs = 1;
            end
            if optimize(1,1) < 0 || optimize(1,1) > 1 || mod(optimize(1,1),1) ~= 0
                fprintf('Fatal Input. Please enter an integer from the list.\n');
                exitFlag = false;
                badexcelinputs = 1;
            end
        end
        switch returnchoice
            % If no coasting is desired for the forward trip,set forward
            % coasting flag to false and all forward coasting parameters
            % to nan because they won't be used.
            case 0
                COAST(1,1) = false;
                departBound(1,1) = nan;
                approachBound(1,1) = nan;
                numIntervals(1,1) = nan;
                coastFraction(1,1) = nan;
                % Allow user to choose forward coasting parameters
            case 1
                COAST(1,1) = true;
                exitFlag = false;
                badexcelinputs = 0;
                if ~optimize(1,1)
                    while ~exitFlag
                        exitFlag = true;
                        if dataEntry == 0
                            departBound(1,1) = input('\nEnter thrust time for departure (as a fraction of one-way trip time): ');
                            approachBound(1,1) = input('Enter thrust time for approach  (as a fraction of one-way trip time): ');
                            coastFraction(1,1) = input('Enter coasting fraction per interval forward: ');
                            numIntervals(1,1) = input('Enter number of coasting intervals forward: ');
                        else
                            if badexcelinputs == 1
                                error('Problem with coasting parameters in excel file.  Please verify that excel all numbers of in the correct range')
                            end
                            departBound(1,1) = xlsread('Input.xlsx','Operation Guidelines','b12:b12');
                            approachBound(1,1) = xlsread('Input.xlsx','Operation Guidelines','b13:b13');
                            coastFraction(1,1) = xlsread('Input.xlsx','Operation Guidelines','b15:b15');
                             numIntervals(1,1) = xlsread('Input.xlsx','Operation Guidelines','b14:b14');
                        end
                        if departBound(1,1) <= 0 || approachBound(1,1) <= 0 || departBound(1,1) >= 1 || approachBound(1,1) >= 1 || coastFraction(1,1) <= 0 || coastFraction(1,1) >= 1 || numIntervals(1,1) <= 0 || mod(numIntervals(1,1),1) ~= 0                            
                            fprintf('\nFatal Inputs. Please enter bounds between 0 and 1\n')
                            exitFlag = false;
                            badexcelinputs = 1;
                        end
                    end
                else
                    numIntervals(1,1) = xlsread('Input.xlsx','Operation Guidelines','b14:b14');
                    departBound(1,1) = nan;
                    approachBound(1,1) = nan;
                    coastFraction(1,1) = nan;
                end
               
        end
end

end