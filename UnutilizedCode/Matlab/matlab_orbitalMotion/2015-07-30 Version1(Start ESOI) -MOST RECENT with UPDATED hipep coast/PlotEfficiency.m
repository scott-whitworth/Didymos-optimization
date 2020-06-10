%PLOTEARTH Plots Earth portion of trip.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PlotAndWriteToFile takes the relevant information from any run of our
%program and plots it. It also saves these plots to the local
%directory with a unique filename. It returns the name of the file
%where it saved the plot from the current run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function PlotEfficiency(thruster,efficiency,inputPower,numberOfEngines)

%Line width for all plots
lineWidth = 3;

%DataPlot is the handle for the figure containing the graphs
DataPlot = figure('Name','EarthPlot','NumberTitle','off','Color','white');
%The plot must be set to be invisible so it does not output a blank graph
%at this point.

sizeArr = 30;
Pin_r = zeros(sizeArr);
eff_r = zeros(sizeArr);

for i = 1:1:sizeArr
    Pin_r(i) = (i+10);
end

% Pin_r = inputPower*numberOfEngines*1000/y(1)^2;
for i = 1:1:sizeArr
    eff_r(i) = efficiency(1)*(Pin_r(i)/(inputPower*numberOfEngines)+efficiency(2)/efficiency(3))^efficiency(4)+efficiency(1)/efficiency(5);
end

subplot(2,1,1)
plot(Pin_r,eff_r,'r','LineWidth',lineWidth)

title(thruster)
xlabel('Power (kW)')
ylabel('Efficiency')

%Gets the next free file name
%File names follow the format 'graph###', where ### is a three-digit number
plotFileName = freename('./','EfficiencyPlot',3);

%Saves file with the generated name
saveas(DataPlot,plotFileName,'fig');

%Closes the figure
close(DataPlot);

end