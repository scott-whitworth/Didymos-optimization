%PLOTEARTH Plots Earth portion of trip.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PlotAndWriteToFile takes the relevant information from any run of our
%program and plots it. It also saves these plots to the local
%directory with a unique filename. It returns the name of the file
%where it saved the plot from the current run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function PlotEarth(order,direction,finalTime,uExhaust,accelInitial,accelFinal,forwardRadius,returnRadius,entryVelocity,entryAngle,asteroidName,thrusterName,powerLevel,thrusterStructs)

maxEscapeTime = 5;

[~,~,~,Tdata,Ydata] = EscapeEarth(order,direction,maxEscapeTime,uExhaust,accelInitial,accelFinal,forwardRadius,returnRadius,entryVelocity,entryAngle);

%Line width for all plots
lineWidth = 3;

%Persistent variables are created to store data on the return trip so that
%both forward and return trips can be plotted when PlotAndWriteToFile is
%called by the forward trip.
persistent TdataReturn;
persistent YdataReturn;
persistent reverseTimeReturn;
persistent finalTimeReturn;

%Converts from non-dimensional time to months;
Tdata = Tdata.*12; 

%Stores time data for return and forward trips seperately, and adds the
%total forward trip time to the return time data so they will plot as a
%continuous trip.
if direction == Direction.FORWARD
    TdataForward = Tdata;
    finalTimeForward = finalTime*12;
    TdataReturn = TdataReturn + finalTimeForward;
else
    TdataReturn = Tdata;
    reverseTimeReturn = TdataReturn;
    finalTimeReturn = finalTime*12;
end

%Convert radial velocity and specific angular momentum to dimensional
%values.
Ydata(:,3) = Ydata(:,3).*(Constants.MCONVERSION/Constants.SCONVERSION);
Ydata(:,4) = Ydata(:,4).*(Constants.MCONVERSION^2/Constants.SCONVERSION);

%DataPlot is the handle for the figure containing the graphs
DataPlot = figure('Name','EarthPlot','NumberTitle','off','Color','white');
%The plot must be set to be invisible so it does not output a blank graph
%at this point.
set(DataPlot,'visible','off');

%Store trip data for the forward and return trip separately.
if direction == Direction.RETURN
    YdataReturn = Ydata;
    YdataForward = [];
else
    YdataForward = Ydata;
end

%Generate plots of thrust angle (alpha), radial velocity (u), and specific 
%angular momentum (h).
if direction == Direction.FORWARD
    set(DataPlot,'visible','on');
    
    %In order to plot return data correctly, the time in which it is
    %plotted with must be reversed.
    TdataReturn = -1.*reverseTimeReturn + (finalTimeReturn+finalTimeForward);
    
    %Plot u, the radial velocity in reference to the sun
    subplot(2,1,1)
    plot(TdataReturn,YdataReturn(:,3),'r','LineWidth',lineWidth)
    hold on
    plot (TdataForward,YdataForward(:,3),'LineWidth',lineWidth)
    title({asteroidName;thrusterName;powerLevel;'Radial Velocity'})
    xlabel('Time(months)')
    ylabel('u(m/s)')
    hold off

    %Plot h, the specific angular momentum in reference to the sun
    subplot(2,1,2)
    plot(TdataReturn,YdataReturn(:,4),'r','LineWidth',lineWidth)
    hold on
    plot(TdataForward,YdataForward(:,4),'LineWidth',lineWidth)
    title('Specific Angular Momentum')
    xlabel('Time(months)')
    ylabel('h(m^2/s)')
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section saves the data plot with a unique filename and then returns
%what that filename is so that it can be accessed later.  It is designed
%this way because this iteration of the program will be iterated through
%hundreds of times overnight.  We want to have all the plots available and
%the reference numbers will be stored along with other data in an Excel
%spreadsheet

%Gets the next free file name
%File names follow the format 'graph###', where ### is a three-digit number
plotFileName = freename('./',thrusterStructs.EarthPlotFileName,3);

if direction == Direction.FORWARD
    %Saves file with the generated name
    saveas(DataPlot,plotFileName,'fig');
    
    %Closes the figure
    close(DataPlot);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create separate figure to hold the plot of the trip trajectory.
ThrustPlot = figure('Name','Trip Trajectory','NumberTitle','off','Color','white');
set(ThrustPlot,'visible','off');
%ArrowPlot() is called to plot the trajectory of the trip, including thrust
%vectors and circles to denote the Earth and asteroid.
PlotTraj(YdataReturn,YdataForward,direction,lineWidth)
ThrustFileName = freename('./',thrusterStructs.EarthTripFileName,3);
if direction == Direction.FORWARD
    set(ThrustPlot,'visible','on');
    %Saves file with the generated name
    saveas(ThrustPlot,ThrustFileName,'fig');
    
    %Close the trip trajectory figure
    close(ThrustPlot);
    
    %Clear persistent variables
    clear TdataReturn;
    clear YdataReturn;
end

end

