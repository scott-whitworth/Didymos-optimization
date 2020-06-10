%INTEGRATESHIPMOTION Integrate ship through heliocentric orbit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculates the initial conditions of the ship in sun sphere
%of influence and then integrates through its heliocentric trip specified
%by time t.  It returns time and position vectors to describe the trip.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by landingObjectiveFunction in fminSearchLanding
function [t,Y] = IntegrateShipMotion(cConverge,t,order,direction,uExhaust,finalMass,efficiency,mdot,numberOfEngines,inputPower,escapeVelocity,escapeEarthAngle,initialConditions,isSolar,object,InitialPositionFlag,timeGuess,interpolateStruct,throttleFlag)

%cConverge(end) is the rotation of earth's reference frame. 0 is defined as
%being where the exit position of the ship is collinear with the radial
%direction of the earth wrt the sun.  See diagram on server.
%cConverge(end) = 0 points away from the sun, and cConverge(end) = pi
%points towards the sun.

if direction == Direction.FORWARD
    %If it is the forward trip and it is at the beginning of the
    %integration it is at the asteroid's location
    if InitialPositionFlag
        y0 = initialConditions;
        y0(5) = 0;
    else
        %If coasting and not at the beginning of firstY the initial
        %conditions are the ship's position from the end of the last
        %integration
        y0 = initialConditions;
    end
    
elseif direction == Direction.RETURN
    if InitialPositionFlag
        %Earth radius + ESOI * cos (rotation of earth frame)
        y0(1) = initialConditions(1)+Constants.ESOI*cos(cConverge(end));
        
        %Initial Angular Position is calculated using law of sines and the
        %triangle made by ESOI Radius, the radius of earth, and the radius from
        %the Sun to the spacecraft. The angle at the suns apex is the
        %calculated and added to the angle of the Earth to calculate the angle
        %of the spacecraft.
        %InitialAngularPosition = AngularPositionEarth + asin(RadialPositionShipEarth * sin((pi - AngleOfRotationOfEarthReferenceFrame) / ShipInitialRadius))        %InitialAngularPosition = AngularPositionEarth + asin(ESOI * sin((Angle in earths frame measure from R_earth)) / ShipInitialRadius wrt Sun))
        y0(2) = initialConditions(2)+asin(sin(pi-cConverge(end))*Constants.ESOI/y0(1));
        
        
        %Calculates initial radial velocity using earths radial velocity, the ship's
        %scalar velocity, and angle of escape relative to Earth.
        %escapeEarthAngle is the escape angle relative to earth and
        %is defined to be 0 when the velocity is entirely angular
        %and 90 when it is entirely radial. This equation calculates the ships
        %initial radius from the sun by combining these values.
        %InitialRadialVelocity = sin(EarthEscapeAngleRelativeToEarth) * EscapeVelocity + RadialVelocityEarth
        y0(3) = sin(escapeEarthAngle)*escapeVelocity + initialConditions(3);
        
        
        %Calculates initial specific angular momementum of ship using earth's
        %specific angular momementum, the ships scalar velocity, escape angle,
        %and initial radius.
        %InitialSpecificAngularMomentum = (cos(EarthEscapeAngleRelativeToEarth) * EscapeVelocity) * ESOI Radius + SpecificAngularMomentumEarth
        y0(4) = (cos(escapeEarthAngle)*escapeVelocity)*Constants.ESOI+Constants.EARTHH;
        
        %Initial mass expended = 0
        y0(5) = 0;
    else
        % The initial conditions are that of the ship.
        y0 = initialConditions;
    end
end
%ode45Fun is used for all orbital integrations in this program. To know if
%its object has thrust and to know what the gravitational center is it
%accesses a global indicator. In this case it is integrating the Ship in
%sun sphere of influence.
% object = Object.SHIP_SUN_SPHERE;

initialStep = timeGuess / 200;

%Set tolerance for heliocentric integration
% options = odeset('RelTol',Constants.ODE45_RELTOL,'AbsTol',Constants.ODE45_ABSTOL);
options = odeset('RelTol',Constants.ODE45_RELTOL,'InitialStep',initialStep,'MaxStep',abs(t(end)-t(1))/2);

%We initialize the struct here in order to save time. if the struct was
%called in SolveOde45 theni would be reconstructed millions of times.
minThrust = interpolateStruct(1).minThrust;
minEfficiency = interpolateStruct(1).minEfficiency;

%Struct created in GetInterpolData
powerInConstMDOT = interpolateStruct(1).powerInConstMDOT;
etaConstMDOT =  interpolateStruct(1).etaConstMDOT;
uExConstMDOT =  interpolateStruct(1).uExConstMDOT;
mDotConstMDOT =  interpolateStruct(1).mDotConstMDOT;

powerInConstUEX = interpolateStruct(1).powerInConstUEX;
etaConstUEX =  interpolateStruct(1).etaConstUEX;
uExConstUEX =  interpolateStruct(1).uExConstUEX;
mDotConstUEX =  interpolateStruct(1).mDotConstUEX;
        
%ode45 gets passed the initial conditions and the time to integrate for and
%returns vectors of time and coordinate values.        
[t,Y] = solveOde45Fun(t,y0,order,direction,object,cConverge,options,uExhaust,NaN,finalMass,efficiency,mdot,numberOfEngines,inputPower,isSolar,throttleFlag,minThrust,minEfficiency,powerInConstMDOT,etaConstMDOT,uExConstMDOT,mDotConstMDOT,powerInConstUEX,etaConstUEX,uExConstUEX,mDotConstUEX);
end
