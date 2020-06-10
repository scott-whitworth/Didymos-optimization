%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function keeps the F vals,Ydata, and Tdata if it is the lowest values so far.  
% This is passed to IntegrateForPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by the landingObjectiveFunction in fminsearch
function [lowestF,MinYdata,MinTdata] = StoreLowestFvalData(F,Ydata,Tdata,lowestF,MinYdata,MinTdata)

%Takes in the current stored lowest data, and then it checks the last Fval
%with the current lowest Fval. If the last Fval is lower it is stored along
%with its corresponding Y and T data

if F < lowestF
    lowestF = F;
    MinYdata = Ydata;
    MinTdata = Tdata;
end
end

