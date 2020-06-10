%CALCULATEENGINEPARAMETERS Calculate global thruster variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is used to calculate accelInitial.  Since the engine is
%always on in our simulation current acceleration can be calculated only
%using accelInitial and current trip time. Most of the global constants are
%initialized in UserInputConstants() or InitializeConstants(). The fuel
%mass for a given trip is passed in in kg and then accelInitial is
%calculated with dimensional calculations, non dimensionalized, and stored
%in a global variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by TimeObjectiveFunction
function [mdot,accelFinal,finalMass] = CalculateEngineParameters(direction,numberOfEngines,uExhaust,inputPower,efficiency,massPowerPlant,massThruster,massPayload,massStruct,massSample)

%All the following calculations are done with dimensional values;

%Power available at LEO for thrusting (kW) per Thruster
Pth = efficiency(1) * inputPower;

%Dimensional mass flow rate (kg/s) per Thruster
mdot = 2*(Pth*1000)/(uExhaust)^2;

%Dimensional thrust (N)
thrust = numberOfEngines * uExhaust * mdot;

%Returns mass of spaceship without fuel (kg) for the return. For the
%forward trip the fuel from the return is included in massPayload. This is
%calculated in MinimizeTripTime
dryMass = massStruct + massPayload + massPowerPlant + massThruster;

if direction == Direction.RETURN
    %Return finalMass is the dryMass and the massSample
    finalMass = dryMass + massSample;
else
    %Forward finalMass is just the dryMass
    finalMass = dryMass;
end

%Dimensional final acceleration at the end of the trip (m/s^2)
accelFinal = thrust/(finalMass);
end
