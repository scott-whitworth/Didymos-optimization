%CALCEARTHCONDITIONS Return final r, theta, and u of Earth.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is used to calculate earths initial conditions. The position
%of the earth at landing is calculated in GetFinalState() which is called
%in InitalizeConstants(). This function accesses those final conditions
%through global variables and takes the trip time as an argument. It then
%reverse integrates the earths motion to find its position at the beginning
%of the trip.  This function is called every time a trip time guess is made
%before the search for an optimum trajectory takes place.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by TimeObjectiveFunction
function [orbitalConditions,t] = CalcOrbitConditions(timeOrbit,order,direction,orbitalElement,uExhaust,closeApp)

%ode45Fun is used to integrate all objects in this program.  To allow
%for that, G_object acts as an indicator to tell the function in what
%sphere of influence the integration is taking place and if there is
%thrust.
object = Object.NO_THRUST;

%Set tolerance integration
% options = odeset('RelTol',Constants.ODE45_RELTOL,'AbsTol',Constants.ODE45_ABSTOL);
options = odeset('RelTol',Constants.ODE45_RELTOL);

%Returns final radial position, angular position, and radial velocity of
%Earth. Since the integration is backwards, the final results are actually
%the initial conditions of our problem.
if orbitalElement == 1
    if direction == Direction.FORWARD
        %Integrates the equations of motion backward in time for the earth
        [~,Y] = solveOde45Fun([timeOrbit,0],closeApp,order,direction,object,0,options,uExhaust,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    elseif direction == Direction.RETURN
        %Integrates the equations of motion forward in time for the earth
        [~,Y] = solveOde45Fun([0,timeOrbit],closeApp,order,direction,object,0,options,uExhaust,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    end
    t = NaN;
    %The end values are returned together in a vector.
    orbitalConditions = [Y(end,1),Y(end,2),Y(end,3),Y(end,4),0];
    
    %Integrates the equations of motion forward in time for the asteroid
elseif orbitalElement == 0
    if direction == Direction.RETURN
        [~,Y] = solveOde45Fun([0,timeOrbit],closeApp,order,direction,object,0,options,uExhaust,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    end
    t = NaN;
    %The end values are returned together in a vector.
    orbitalConditions = [Y(end,1),Y(end,2),Y(end,3),Y(end,4),0];
    
    %Integrates the equations of motion forward in time during the spaceship's
    %waiting peroid
elseif orbitalElement == 2
    if direction == Direction.FORWARD
        [t,Y] = solveOde45Fun([0,timeOrbit],closeApp,order,direction,object,0,options,uExhaust,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    end
    %All the values are returned together in a vector.
    orbitalConditions = Y;
end

end
