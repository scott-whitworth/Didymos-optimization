%PLOTANDWRITETOFILE Plots alpha, u, and h. Calls ArrowPlot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PlotAndWriteToFile takes the relevant information from any run of our
%program and plots it. It also saves these plots to the local
%directory with a unique filename. It returns the name of the file
%where it saved the plot from the current run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function [plotFileName,rError,thetaError,uError,hError] = PlotAndWriteToFile(order,direction,phiCoeff,fval,escapeTime,finalTime,timeWait,uExhaust,asteroidConditionsForward,asteroidConditionsReturn,asteroidName,powerLevel,thrusterName,Tdata,Ydata,thrusterStructs,Adata,departData,coastData,thrustData,approachData,COAST,LegNumber,AsteroidSemimajor,AsteroidEccentricity,MdotData,TmdotData)

%Line width for all plots
lineWidth = 3;

%Persistent variables are created to store data on the return trip so that
%both forward and return trips can be plotted when PlotAndWriteToFile is
%called by the forward trip.
persistent TdataReturn;
persistent YdataReturn;
persistent AdataReturn;
persistent MdotdataReturn;
persistent TmdotdataReturn;
persistent AdataReturnFlipped;
persistent YdataDepartReturn;
persistent YdataApproachReturn;
persistent YdataCoastReturn;
persistent YdataThrustReturn;
persistent alphaReturn;
persistent gammaReturn;
persistent finalTimeReturn;
persistent returnEscapeTime;
persistent AdataForward;
persistent MdotdataForward;
persistent TmdotdataForward;
%Make variable names more readable
if direction == Direction.FORWARD
    goalRadius = asteroidConditionsForward(1);
    goalTheta = asteroidConditionsForward(2);
    goalU = asteroidConditionsForward(3);
    goalH = asteroidConditionsForward(4);
elseif direction == Direction.RETURN
    goalRadius = asteroidConditionsReturn(1);
    goalTheta = asteroidConditionsReturn(2);
    goalU = asteroidConditionsReturn(3);
    goalH = asteroidConditionsReturn(4);
end

%Integrate the ships motion while it is coasting
if direction == Direction.FORWARD
    waitData = zeros(5,1);
    waitData(1) = Ydata(1,1);
    waitData(2) = Ydata(1,2);
    waitData(3) = Ydata(1,3);
    waitData(4) = Ydata(1,4);
    
    orbitalElement = 2;
    if timeWait > 0
        [YdataWait,TdataWait] = CalcOrbitConditions(timeWait,order,direction,orbitalElement,uExhaust,waitData);
        TdataWait = (TdataWait + finalTime).*12;
        timeWait = timeWait*12;
    else
        YdataWait = zeros(4,1);
    end
elseif direction == Direction.RETURN
    %Initialize to prevent undefined variable errors
    TdataWait = 0;
    YdataWait = 0;
end

finalRadius = Ydata(end,1);
rError = abs(finalRadius-goalRadius);
finalTheta = Ydata(end,2);
thetaError = abs(finalTheta-goalTheta);
finalU = Ydata(end,3);
uError = abs(finalU-goalU);
finalH = Ydata(end,4);
hError = abs(finalH-goalH);

%The following displays the dimensionalized conditions of the ship and
%asteroid at the time of landing/takeoff. It calculates the differences
%between these values and displays them as well. This is done to see how
%closely the objective function is met.
Rocket_rad = finalRadius*Constants.MCONVERSION;
Asteroid_rad = goalRadius*Constants.MCONVERSION;
Dif_rad = abs(Asteroid_rad - Rocket_rad);
Dif_circ = thetaError*Rocket_rad;

%Convert to (m/s);
Rocket_U = finalU*(Constants.MCONVERSION/Constants.SCONVERSION);
Asteroid_U = goalU*(Constants.MCONVERSION/Constants.SCONVERSION);
Dif_U = abs(Asteroid_U - Rocket_U);

%Convert to (m^2/s)
Rocket_h = finalH*(Constants.MCONVERSION^2/Constants.SCONVERSION);
Asteroid_h = goalH*(Constants.MCONVERSION^2/Constants.SCONVERSION);
Dif_h = abs(Asteroid_h - Rocket_h);
Dif_v_theta = Dif_h/Rocket_rad;

%We don't pass the entire phiCoeff vector because the last element
%represents the escape angle from earth and should not be evaluated as part
%of the fourier series.
PhiData = EvaluateFourierSeries(direction,phiCoeff(1:(end-1)),Tdata);

%Set betaData, the angle between the heliocentric tangent and the velocity
%vector of the spacecraft.

betaData = atan((Ydata(:,3).*Ydata(:,1).^2)./Ydata(:,4));
betaData = betaData.';

%Set alphaData, the angle between the velocity vector of the spacecraft and
%the thrust vector of the spacecraft.
alphaData = mod(PhiData - betaData,2*pi);
for i = 1:size(alphaData,2)
    if (alphaData(i) > pi)
        alphaData(i) = mod(alphaData(i),pi) - pi;
    end
end

%Set gammaData. At each point along the trajectory of the spacecraft, if a
%cartesian coordinate system is created with its origin at the current
%point along the trajectory, then gamma is the angle between the positive
%x-axis and thrust vector.  This allows for plotting of the thrust vectors.
gammaData = (pi/2) + Ydata(:,2).' - PhiData;

%Converts from non-dimensional time to months;
Tdata = Tdata.*12;
TmdotData = TmdotData.*12;
escapeTime = escapeTime.*12;

%Stores time data for return and forward trips seperately, and adds the
%total forward trip time to the return time data so they will plot as a
%continuous trip.
if direction == Direction.FORWARD
    TdataForward = Tdata;
    finalTimeForward = finalTime.*12;
    TdataReturn = TdataReturn + finalTimeForward + timeWait;
    TmdotdataReturn = TmdotdataReturn + finalTimeForward + timeWait;
else
    returnEscapeTime = escapeTime;
    TdataReturn = Tdata;
    finalTimeReturn = finalTime.*12;
end

%Convert radial velocity and specific angular momentum to dimensional
%values.
Ydata(:,3) = Ydata(:,3).*(Constants.MCONVERSION/Constants.SCONVERSION);
Ydata(:,4) = Ydata(:,4).*(Constants.MCONVERSION^2/Constants.SCONVERSION);

%Various values from the trip data are converted to string to be inserted
%into the title of the plot.
uexname = num2str(uExhaust);
fvalname = num2str(fval);
rgoalname = num2str(goalRadius);
ordername = num2str((length(phiCoeff)-2)/2);
timename = num2str(finalTime);
escapeanglename = num2str(phiCoeff(end));
name = [' fval ',fvalname,' rf ',rgoalname,' order ',ordername,' time ',timename,' escape angle ',escapeanglename];

%DataPlot is the handle for the figure containing the graphs
DataPlot = figure('Name',name,'NumberTitle','off','Color','white');
%The plot must be set to be invisible so it does not output a blank graph
%at this point.
set(DataPlot,'visible','off');

%Store trip data for the forward and return trip separately.
if direction == Direction.RETURN
    YdataReturn = Ydata;
    AdataReturn = Adata;
    AdataReturnFlipped = AdataReturn;
    YdataDepartReturn = departData;
    YdataApproachReturn = approachData;
    YdataCoastReturn = coastData;
    YdataThrustReturn = thrustData;
    MdotdataReturn = MdotData(:,1)./Constants.SCONVERSION;
    TmdotdataReturn = TmdotData;
    alphaReturn = alphaData;
    gammaReturn = gammaData;
    YdataForward = [];
    YdataDepartForward = [];
    YdataApproachForward = [];
    YdataCoastForward = [];
    YdataThrustForward = [];
    alphaForward = [];
    gammaForward = [];
    AdataForward = [];
else
    YdataForward = Ydata;
    AdataForward = Adata;
    YdataDepartForward = departData;
    YdataApproachForward = approachData;
    YdataCoastForward = coastData;
    YdataThrustForward = thrustData;
    MdotdataForward = MdotData(:,1)./Constants.SCONVERSION;
    TmdotdataForward = TmdotData;
    alphaForward = alphaData;
    gammaForward = gammaData;
end

%Generate plots of thrust angle (alpha), radial velocity (u), and specific
%angular momentum (h).
if direction == Direction.FORWARD
    waitThrust = 1:.1:timeWait;
    waitThrustTime = waitThrust;
    waitThrust = zeros(size(waitThrust));
    earthThrust = 1:.1:escapeTime*12;
    earthThrustTime = earthThrust;
    earthThrust = zeros(size(earthThrust));
    
    set(DataPlot,'visible','on');
    
    %Plot u, the radial velocity in reference to the sun
    subplot(2,2,1)
    plot(TdataReturn,YdataReturn(:,3),'r.')
    hold on
    plot(TdataForward,YdataForward(:,3),'.')
    axis([0 (finalTimeReturn+returnEscapeTime+timeWait+finalTimeForward) -Inf Inf])
    title('Radial Velocity')
    xlabel('Time(months)')
    ylabel('u(m/s)')
    hold off
    
    %Plot h, the specific angular momentum in reference to the sun
    subplot(2,2,3)
    plot(TdataReturn,YdataReturn(:,4),'r.')
    hold on
    plot(TdataForward,YdataForward(:,4),'.')
    axis([0 (finalTimeReturn+returnEscapeTime+timeWait+finalTimeForward) -Inf Inf])
    title('Specific Angular Momentum')
    xlabel('Time(months)')
    ylabel('h(m^2/s)')
    hold off
    
    if ~COAST(1,1) 
        AdataReturnPlot = AdataReturnFlipped;
    else
        AdataReturnPlot = AdataReturn;
    end
    
    %Plot acceleration due to thrusting
    subplot(2,2,2)
    plot(TdataReturn,AdataReturnPlot,'r.')
    hold on
    plot(TdataForward,AdataForward,'b.')
    axis([0 (finalTimeReturn+returnEscapeTime+timeWait+finalTimeForward) -Inf Inf])
    title('Thrust Acceleration')
    xlabel('Time(months)')
    ylabel('a(m/s^2)')
    hold off

    %Plot mdot, mass flow rate
    subplot(2,2,4)
    plot(TmdotdataReturn,MdotdataReturn,'r.')
    hold on
    plot(TmdotdataForward,MdotdataForward,'.')
    axis([0 (finalTimeReturn+returnEscapeTime+timeWait+finalTimeForward) 0 Inf])
    title('Propellant Flow Rate')
    xlabel('Time(months)')
    ylabel('mdot(kg/s)')
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
plotFileName = freename('./',thrusterStructs.HelioPlotFileName,3);

if direction == Direction.FORWARD
    %Saves file with the generated name
    saveas(DataPlot,plotFileName,'fig');
    
    %Closes the figure
    close(DataPlot);
end

if direction == Direction.FORWARD
    % Plots the phi coefficients in its own file
    
    PhiPlot = figure('Name',name,'NumberTitle','off','Color','white');
    set(PhiPlot,'visible','on');
    %Plot alpha, the thrust angle in reference to the space craft
    plot(TdataReturn,alphaReturn,'r.')
    hold on
    plot(TdataForward,alphaForward,'.')
    %Set y-axis to plot from -pi to pi.
    axis([0 (finalTimeReturn+returnEscapeTime+timeWait+finalTimeForward) -pi pi])
    set(gca,'YTick',-pi:pi/2:pi)
    set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
    title({asteroidName;thrusterName;powerLevel;'Thrust Angle (alpha)'})
    xlabel('Time(months)')
    ylabel('Thrust Angle (rad)')
    hold off
    
    fileName = strcat(thrusterStructs.HelioPlotFileName,'PhiPlot');
    PhiplotFileName = freename('./',fileName,3);
    %Saves file with the generated name
    saveas(PhiPlot,PhiplotFileName,'fig');
    
    %Closes the figure
    close(PhiPlot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create separate figure to hold the plot of the trip trajectory.
ThrustPlot = figure('Name','Trip Trajectory','NumberTitle','off','Color','white');
set(ThrustPlot,'visible','off');
%ArrowPlot() is called to plot the trajectory of the trip, including thrust
%vectors and circles to denote the Earth and asteroid.

TripTrajectoryPlot(YdataReturn,YdataForward,YdataDepartReturn,YdataApproachReturn,YdataCoastReturn,YdataThrustReturn,YdataDepartForward,YdataApproachForward,YdataCoastForward,YdataThrustForward,YdataWait,timeWait,gammaReturn,gammaForward,direction,lineWidth,asteroidName,powerLevel,thrusterName,AdataReturn,AdataForward,COAST,AsteroidSemimajor,AsteroidEccentricity)

if direction == Direction.FORWARD
    set(ThrustPlot,'visible','on');
    %Saves file with the generated name
    ThrustFileName = freename('./',thrusterStructs.HelioTripFileName,3);
    saveas(ThrustPlot,ThrustFileName,'fig');
    
    %Close the trip trajectory figure
    close(ThrustPlot);
    
    %Clear persistent variables
    clear phiReturn;
    clear TdataReturn;
    clear YdataReturn;
    clear AdataReturn;
    clear AdataReturnFlipped;
    clear YdataDepartReturn;
    clear YdataApproachReturn;
    clear YdataCoastReturn;
    clear YdataThrustReturn;
    clear alphaReturn;
    clear gammaReturn;
    clear finalTimeReturn;
    clear returnEscapeTime;
    clear AdataForward;   

end
end
