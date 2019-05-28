%PLOTTRAJ Plots the path of the space craft.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create persistent variables to store plot data on the return trip in order
%to be plotted when ArrowPlot is called by the forward trip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by PlotAndWriteToFile
function PlotTraj(YdataReturn,YdataForward,direction,lineWidth)

%Persistent variables are created to store data on the return trip so that
%both forward and return trips can be plotted when PlotAndWriteToFile is
%called by the forward trip.
persistent xReturn;
persistent yReturn;

%Set properties of plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Size of circles plotted
circleSize = 50;
%Color of circle denoting Earth
earthGreen = [0 .5 0];
%Color of circle denoting the asteroid
asteroidOrange = [.8 .4 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the trajectory of the space craft
if direction == Direction.FORWARD
    %Plot forward trajectory
    forwardPlot = polar(YdataForward(:,2),YdataForward(:,1));
    %Plot is held open so subsequent calls to plot function will plot on
    %the same figure.
    hold on
    %Plot return trajectory
    returnPlot = polar(YdataReturn(:,2),YdataReturn(:,1),'r');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %removes angle labels from the polar plot
    hHiddenText = findall(gca,'type','text');
    k = 1;
    hObjToDelete = zeros(12,1);
    for ang = 0:30:330
        hObjToDelete(k) = findall(hHiddenText,'string',num2str(ang));
        k = k + 1;
    end
    delete(hObjToDelete);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Polar plots cannot take standard LineSpec arguments, so the line width
    %must be set using set(plotHandle,arguments...).
    set(returnPlot,'LineWidth',lineWidth)
    set(forwardPlot,'LineWidth',lineWidth)
    %Set title of plot
    title({'Plot of trajectory with vectors denoting thrust angle';'(Circles denote radius from Earth in AU)'})
end

%The trajectory data for the return and forward trips must be converted to
%cartesian coordinates for graphing the thrust vectors
if direction == Direction.RETURN
    %Convert the points in the trajectory from polar to cartesian
    %coordinates. x and y are used as the points of origin for the thrust
    %vectors.
    [xReturn,yReturn] = pol2cart(YdataReturn(:,2),YdataReturn(:,1));
    %The x and y arrays must be transposed in order to match the size of u
    %and v
    xReturn = xReturn.';
    yReturn = yReturn.';
else
    %Convert the points in the trajectory from polar to cartesian
    %coordinates. x and y are used as the points of origin for the thrust
    %vectors.
    [xForward,yForward] = pol2cart(YdataForward(:,2),YdataForward(:,1));
    %The x and y arrays must be transposed in order to match the size of u
    %and v
    xForward = xForward.';
    yForward = yForward.';
end

%Plot circles marking the Earth, in green, and the asteroid, in orange.
if direction == Direction.FORWARD
    %Plot orange circle at the start of the return trajectory
    scatter(xReturn(end),yReturn(end),circleSize,asteroidOrange,'LineWidth',lineWidth);
    %Plot green circle at the end of the return trajectory
    scatter(xReturn(1),yReturn(1),circleSize,earthGreen,'LineWidth',lineWidth);
    %Plot green circle at the start of the forward trajectory
    scatter(xForward(1),yForward(1),circleSize,earthGreen,'LineWidth',lineWidth);
    %Plot orange circle at the end of the forward trajectory
    scatter(xForward(end),yForward(end),circleSize,asteroidOrange,'LineWidth',lineWidth);
    
    %Persistent variables are cleared
    clear xReturn;
    clear yReturn;
end

%Close plot
hold off
end