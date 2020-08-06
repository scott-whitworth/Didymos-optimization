%This function is designed to take the new Pin values that were created by
%GetInterpolData and plot it

function EffPlot( filename, chooseExtrap, extrapPercent, mdot )

DataPlot = figure('Name','EfficiencyPlot','Numbertitle','off','Color','white');
%makes the plot invisible until the plot is writen

%Calls GetInterpolData to have the PowerVector in the function. The
%PowerVector will be called later so that we can calculate the efficiency.
[PowerVector, uExhaustVector] = GetInterpolData(filename,chooseExtrap,extrapPercent)

%Sizes both the PowerVector and the uExhaustVector matricies in order to
%run a for loop. At the Moment the for loop is commented out becasue it did
%not do what I wanted it to.
uExSize = length(uExhaustVector);
powerSize = length(PowerVector);

%powerUEx = [PowerVector,uExhaustVector];

eta = zeros(powerSize,1);

for i = 1:1:powerSize-1
    eta(i,1) = ((.5*mdot*((uExhaustVector(i,1))^2))/(PowerVector(i,1)));
end

%Plots the functions with the PowerVector in the X and eta in the Y
plot(PowerVector,eta,'g.')
%plot(PowerVector,uExhaustVector,'g.')

%Gives the Plot that the code made a title and lables
title(filename)
xlabel('Power')
ylabel('Efficiency')

%plotFileName will name the plot for the code with a 3 digit end so that it
%doesn't write over other plots.
plotFileName = freename('./','EfficiencyPlot',3);

%saves Plot
saveas(DataPlot,plotFileName,'fig');

%close(DataPlot);
eta

end

