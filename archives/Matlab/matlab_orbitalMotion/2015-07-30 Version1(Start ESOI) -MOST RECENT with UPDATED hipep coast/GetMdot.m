%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function calculates the value of total mdot at different locations in the
%journey based on the slope of y5. This is calculated for plotting
%purposes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MdotData,TmdotData]= GetMdot(Mdata,Tdata,numberOfEngines,mdot)

%Calculating slope so number of points in 1 less then in Mdata
numPoints = (size(Mdata)-1);

%Initialize vectors
MdotData = zeros(numPoints(1,1),1);
TmdotData = zeros(numPoints(1,1),1);

for i=1:1:numPoints(1,1)
    
    %Calculate slope of y5
    PreMdotData(i,1) = (Mdata(i+1,1) - Mdata(i,1))/(Tdata(i+1,1)-Tdata(i,1));
    
    %Calculate associated time values
    TmdotData(i,1) = (Tdata(i+1,1) + Tdata(i,1))/2;
    
    
    MdotData(i,1) = PreMdotData(i,1);
    %This cap is for if the slope provides an mdot greater than the maximum
    %mdot for the engine it sets its value at that maximum.
    
    %If mdot calculated is greater then maximum mdot, set it to the maximum
%     %mdot value
%     if (PreMdotData(i,1)/Constants.SCONVERSION) > (mdot * numberOfEngines)
%         MdotData(i,1) = mdot * numberOfEngines * Constants.SCONVERSION;
%     else
%         MdotData(i,1) = PreMdotData(i,1);
%     end
end
end
