function [flagStore] = ChooseMinimize()
fprintf('\nOptimize Time or Fuel?\n')
fprintf('0: Time\n')
fprintf('1: Fuel\n')
choice = input(': ');
switch choice
    case 0
        flagStore = true;
    case 1
        flagStore = false;
end
    
end