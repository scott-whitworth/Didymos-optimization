%ESCAPEEARTH Returns the time to escape Earth's SOI for a given thruster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EscapeEarth is passed a time that is guaranteed to be greater than the
%time required to escape earth.  It then starts the ship in low earth orbit
%and has it thrust parallel to its motion.  The integration has a
%termination condition when the ship is outside earth sphere of influence.
%The time at which it escaped earth influence is returned and the angle and
%velocity at which it escaped are stored in global variables for later use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by TimeObjectiveFunction
function [escapeTime,escapeVelocity,escapeEarthAngle,t,Y] = EscapeEarth(order,direction,maxEscapeTime,uExhaust,accelFinal,forwardRadius,returnRadius,entryVelocity,entryAngle,mdot,numberOfEngines)

object = Object.SHIP_EARTH_SPHERE;

%Sets the event condition which stops integration when radius is equal to
%ESOI.
% options = odeset('Events', @EscapeEvent,'RelTol',Constants.ODE45_RELTOL,'AbsTol',Constants.ODE45_ABSTOL);
options = odeset('Events', @EscapeEvent,'RelTol',Constants.ODE45_RELTOL);

%Initializes the vector of initial conditions
yInitial = zeros(5,1);

%The return to earth is modeled by attempting a direct entry into earth's
%atmosphere, whereas the forward escape earth event is a gradual spiral.
if direction == Direction.RETURN
    %Final r (non-dimensional)
    yInitial(1) = returnRadius;
    %Final theta (rad)
    yInitial(2) = 0;
    
    %The two equations
    %entryAngle = atan(u*r/h)
    %entryVelocity = sqrt(u^2+(h/r)^2)
    %were manipulated and solved for u and h in order to find the final
    %values of u and h.
    
    %Final u (non-dimensional)
    yInitial(3) = entryVelocity / sqrt(1 + cot(entryAngle)^2);
    %Final h (non-dimensional)
    yInitial(4) = (yInitial(3) * yInitial(1)) / tan(entryAngle);
    yInitial(5) = 0;
else
    %Initial r (non-dimensional)
    yInitial(1) = forwardRadius;
    %Initial theta (rad)
    yInitial(2) = 0;
    %Initial u (non-dimensional)
    yInitial(3) = 0;
    %Initial h (non-dimensional)
    yInitial(4) = sqrt(Constants.MUE*forwardRadius);
    yInitial(5) = 0;
end

%Integrate equations of motion for Earth escape
[t,Y,~,~,~] = solveOde45Fun([0,maxEscapeTime],yInitial,order,direction,object,0,options,uExhaust,accelFinal,NaN,NaN,mdot,numberOfEngines,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);

%Use pythagorean theorem and the radial and angular velocities
escapeVelocity = sqrt(Y(end,3)^2+(Y(end,4)/Y(end,1))^2);

%atan of (radialVelocity/angularVelocity) defines 0 to be completely
%angular and 90 degrees to be completely radial.
escapeEarthAngle = atan((Y(end,1)*Y(end,3)/Y(end,4)));

%The integration terminates when earth is escaped so the last value will be
%the escape time
escapeTime = t(end);

%C3 gives the total energy of the space craft at escape. This should be
%positive if the space craft actually has enough energy to escape
%Earth. Also, it should be close to zero in order to not expend more
%energy than is necessary. As a non-dimensional number it gives a rough
%percentage of how much extra energy the space craft has at Earth
%escape.
C3 = (1/2)*(Y(end,3)^2 + (Y(end,4)/Y(end,1))^2) - (Constants.G*Constants.ME)/Y(end,1);

%If radial position never reaches escape radius then the function has
%failed and returning a non-number value is the easiest way to demonstrate
%that there was an error.
if Y(end,1)<Constants.ESOI
    fprintf ('EscapeEarth() NEVER REACHED GOAL RADIUS, CHANGE THE VALUE WHICH IT IS PASSED IN THE CODE TO BE GREATER\n')
    escapeTime = nan;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Contains the event which triggers upon switching spheres of influence
%**************************************************************
%TODO
%Something we still need to do is change this so that ESOI is not a
%constant and instead varies with the current distance from the earth to
%the sun. Does this function have to re-integrate earths motion for every
%function call?  That could get very expensive.
%******************************************************************
%Function which handles event.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = EscapeEvent(~,y)

%Looks for current radius to be equal to the point where we switch spheres
%of influence
value(1) = (y(1)-Constants.ESOI);
%Tells MATLAB to stop integrating
isterminal(1) = 1;
%Doesn't care which direction the function is moving
direction(1) = 0;
end
