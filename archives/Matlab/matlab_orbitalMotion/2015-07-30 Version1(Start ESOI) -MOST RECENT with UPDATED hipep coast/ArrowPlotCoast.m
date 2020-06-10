%ARROWPLOT Plots the path of the space craft with vectors to denote thrust.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create persistent variables to store plot data on the return trip in order
%to be plotted when ArrowPlot is called by the forward trip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by PlotAndWriteToFile
function [sizeArrReturn,xReturn,yReturn,uReturn,vReturn] = ArrowPlotCoast(YdataReturn,YdataForward,YdataDepartReturn,YdataApproachReturn,YdataCoastReturn,YdataThrustReturn,YdataDepartForward,YdataApproachForward,YdataCoastForward,YdataThrustForward,YdataWait,timeWait,gammaReturn,gammaForward,direction,lineWidth,asteroidName,powerLevel,thrusterName,AdataReturn,AdataForward,sizeArrReturn,xReturn,yReturn,uReturn,vReturn)


%Set properties of plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Length of thrust vector arrows
arrowSize = 0.1;
%Weight of arrow line
arrowLineWidth = 2;
%Size of circles plotted
circleSize = 50;
%Color of circle denoting Earth
earthGreen = [0 .5 0];
%Color of circle denoting the asteroid
asteroidOrange = [.8 .4 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the trajectory of the space craft
if direction == Direction.FORWARD
    %Plot return trajectory
    returnDepartPlot = polar(YdataDepartReturn(:,2),YdataDepartReturn(:,1),'r');
    hold on
    returnCoastPlot =  polar(YdataCoastReturn(:,2),YdataCoastReturn(:,1),'c.');
    hold on
    returnThrustPlot =  polar(YdataThrustReturn(:,2),YdataThrustReturn(:,1),'r.');
    hold on
    returnApproachPlot =  polar(YdataApproachReturn(:,2),YdataApproachReturn(:,1),'r');
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
    %Plot is held open so subsequent calls to plot function will plot on
    %the same figure.
    hold on
    %Plot forward trajectory
    forwardDepartPlot = polar(YdataDepartForward(:,2),YdataDepartForward(:,1));
    hold on
    forwardCoastPlot =  polar(YdataCoastForward(:,2),YdataCoastForward(:,1),'c.');
    hold on
    forwardThrustPlot =  polar(YdataThrustForward(:,2),YdataThrustForward(:,1),'.');
    hold on
    forwardApproachPlot =  polar(YdataApproachForward(:,2),YdataApproachForward(:,1));
    %Plot waiting time on asteroid
    if timeWait > 0
        waitPlot = polar(YdataWait(:,2),YdataWait(:,1),'k');
    end
    %Polar plots cannot take standard LineSpec arguments, so the line width
    %must be set using set(plotHandle,arguments...).
    set(returnDepartPlot,'LineWidth',lineWidth)
    set(forwardDepartPlot,'LineWidth',lineWidth)
    set(returnCoastPlot,'LineWidth',lineWidth)
    set(forwardCoastPlot,'LineWidth',lineWidth)
    set(returnThrustPlot,'LineWidth',lineWidth)
    set(forwardThrustPlot,'LineWidth',lineWidth)
    set(returnApproachPlot,'LineWidth',lineWidth)
    set(forwardApproachPlot,'LineWidth',lineWidth)
    if timeWait > 0
        set(waitPlot,'LineWidth',lineWidth)
    end
    %Set title of plot
    title({asteroidName;thrusterName;powerLevel;'Plot of trajectory with vectors denoting thrust angle';'(Circles denote radius from sun in AU)'})
end

%The trajectory data for the return and forward trips must be converted to
%cartesian coordinates for graphing the thrust vectors
if direction == Direction.RETURN
    %Convert the points in the trajectory from polar to cartesian
    %coordinates. x and y are used as the points of origin for the thrust
    %vectors.
    [xReturn,yReturn] = pol2cart(YdataReturn(:,2),YdataReturn(:,1));
    %Convert the thrust angle for each point in the trajectory from polar to
    %cartesian coordinates. u and v are thus the components of the thrust
    %vector.
    [uReturn,vReturn] = pol2cart(gammaReturn,repmat(arrowSize,size(gammaReturn,1),1));
    %The x and y arrays must be transposed in order to match the size of u
    %and v
    xReturn = xReturn.';
    yReturn = yReturn.';
    %The variable sizeArr... is defined as the length of x, which is equal to the
    %lengths of y, u, and v
    sizeArrReturn = size(xReturn,2);
else
    %Convert the points in the trajectory from polar to cartesian
    %coordinates. x and y are used as the points of origin for the thrust
    %vectors.
    [xForward,yForward] = pol2cart(YdataForward(:,2),YdataForward(:,1));
    %Convert the thrust angle for each point in the trajectory from polar to
    %cartesian coordinates. u and v are thus the components of the thrust
    %vector.
    [uForward,vForward] = pol2cart(gammaForward,repmat(arrowSize,size(gammaForward,1),1));
    %The x and y arrays must be transposed in order to match the size of u
    %and v
    xForward = xForward.';
    yForward = yForward.';
    %The variable sizeArr... is defined as the length of x, which is equal to the
    %lengths of y, u, and v
    sizeArrForward = size(xForward,2);
end

%Plot circles marking the Earth, in green, and the asteroid, in orange.
if direction == Direction.FORWARD
    %Plot green circle at the end of the return trajectory
    scatter(xReturn(1),yReturn(1),circleSize,earthGreen,'LineWidth',lineWidth);
    %Plot orange circle at the start of the return trajectory
    scatter(xReturn(end),yReturn(end),circleSize,asteroidOrange,'LineWidth',lineWidth);
    %Plot green circle at the start of the forward trajectory
    scatter(xForward(1),yForward(1),circleSize,earthGreen,'LineWidth',lineWidth);
    %Plot orange circle at the end of the forward trajectory
    scatter(xForward(end),yForward(end),circleSize,asteroidOrange,'LineWidth',lineWidth);
end

%The output vectors are thinned to only display at intervals greater than
%the set angle step.
%This is accomplished by placing the indices of the vectors which are
%spaced at a distance less than radStep into the array ind.  These
%indices are then removed from x, y, u, and v so that they will not be
%plotted.
%This must be done for both the return and forward trip data.
radStep = 0.1;
lastVector = 1;
indCounter = 1;
if direction == Direction.FORWARD
    %Forward trip vectors are thinned.
    %Data which is spaced too closely is found and their indices are
    %stored in the array ind
    for i = 2:sizeArrForward
        if (abs(YdataForward(i,2)-YdataForward(lastVector,2)) <= radStep)
            ind(indCounter) = i;
            indCounter = indCounter + 1;
        elseif (abs(AdataForward(i,1)) == 0)
            ind(indCounter) = i;
            indCounter = indCounter + 1;
        else
            lastVector = i;
        end
    end
    %The data stored at the indices in the array ind is removed from x, y,
    %u, and v
    xForward(ind) = [];
    yForward(ind) = [];
    uForward(ind) = [];
    vForward(ind) = [];
    
    %The thrust vectors are plotted as arrows along the trajectory of the
    %space craft by the function quiver(x,y,u,v,...).  The variables x and
    %y denote the location of the origin of the vector while u and v denote
    %the components of the vector. Basic Matlab LineSpec properties also
    %apply.
    quiver(xForward,yForward,uForward,vForward,'k','AutoScale','off','LineWidth',arrowLineWidth)
    quiver(xReturn,yReturn,uReturn,vReturn,'k','AutoScale','off','LineWidth',arrowLineWidth)   
    
else
    %Return trip vectors are thinned.
    %Data which is spaced too closely is found and their indices are
    %stored in the array ind
    for i = 2:sizeArrReturn
        if (abs(YdataReturn(i,2)-YdataReturn(lastVector,2)) <= radStep)
            ind(indCounter) = i;
            indCounter = indCounter + 1;
        elseif (abs(AdataReturn(i,1)) == 0)
            ind(indCounter) = i;
            indCounter = indCounter + 1;
        else
            lastVector = i;
        end
    end
    %The data stored at the indices in the array ind is removed from x, y,
    %u, and v
    xReturn(ind) = [];
    yReturn(ind) = [];
    uReturn(ind) = [];
    vReturn(ind) = [];
 
end
%Close plot
hold off
end