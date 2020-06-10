% Allow user to choose the Payload mass from a preset list.

function Payload = ChoosePayload()

exitFlag = false;
while ~exitFlag
    exitFlag = true;
    fprintf('\nSelect Payload mass (kg)\n')
    fprintf('1  625\n')
    fprintf('2  750\n')
    fprintf('3  875\n\n')
    
    choice = input(': ');
    if choice < 0 || choice > 3 || mod(choice,1) ~= 0
        fprintf('Fatal Input. Please enter an integer from the list.\n');
        exitFlag = false;
    end
end
switch choice
    case 1
        Payload = 625;
    case 2
        Payload = 750;
    case 3
        Payload = 875;

end

end