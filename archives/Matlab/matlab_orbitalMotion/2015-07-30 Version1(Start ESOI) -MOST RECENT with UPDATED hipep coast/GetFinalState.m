%GETFINALSTATE Get position/velocity vectors of Earth/Asteroid at landing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculates the positions of the asteroid and earth at
%landing. It calls GetCloseApprach to recieve their x,y,dx,dy values at
%that time and then converts those to r,u, and theta values.  It also
%calculates H based off of C and the mu of the Sun. The H of earth should
%be a constant throughout all runs of the program but we thought we should
%show how value was calculated and for consistency/symmetry sake we put it
%here with the other calculations for final values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function [earthCloseApp,asteroidCloseApp,AsteroidSemimajor,AsteroidEccentricity] = GetFinalState(astDes)

%Gets vectors for earth and asteroid relative positions at the point of
%closest approach.  The vectors contain x,y,dx,dy in that order.
[AData,EData,AsteroidEccentricity,AsteroidSemimajor] = GetCloseApproach(astDes);

%Eccentricity and semimajor are non-dimentionalized in the text file
asteroidC = AsteroidSemimajor*(1-AsteroidEccentricity^2);

%H represents angular momentum per mass
%H value for asteroid (non-dimensional)
AsteroidH = sqrt(asteroidC*Constants.MUS);    %H value for asteroid (non-dimensional)

%The program can easily adopt any reference system in terms of theta. It
%only requires a final theta value for the earth and asteroid, where
%theta=0 is defined is irrelevant. For simplicities sake we will put
%theta=0 on the x axis of whatever coordinate grid corresponds with the x
%and y values being read in from the text file. If it ever becomes
%important to adopt a different reference system it should be easy to
%change.

%ASR = Asteroid-Sun Radius
[ATheta,ASR] = cart2pol(AData(1),AData(2)); %Converts x,y values to r, theta with theta=0 defined on x axis
%ESR = Earth-Sun Radius
[ETheta,ESR] = cart2pol(EData(1),EData(2)); %Converts x,y values to r, theta with theta=0 defined on x axis

%Object radial velocity.
%The ratio between the x and y coordinates and
%the radius is used to weight the velocity in x and y, the two are then
%added to get radial velocity. The equation can be derived by taking the
%derivative of r=(x^2+y^2)^.5
earthU = EData(3)*EData(1)/ESR+EData(4)*EData(2)/ESR;
asteroidU = AData(3)*AData(1)/ASR+AData(4)*AData(2)/ASR;

%Return conditions in arrays
earthCloseApp = zeros(1,5);
earthCloseApp(1,1) = ESR;
earthCloseApp(1,2) = ETheta;
earthCloseApp(1,3) = earthU;
earthCloseApp(1,4) = Constants.EARTHH;
earthCloseApp(1,5) = 0; 

asteroidCloseApp = zeros(1,5);
asteroidCloseApp(1,1) = ASR;
asteroidCloseApp(1,2) = ATheta;
asteroidCloseApp(1,3) = asteroidU;
asteroidCloseApp(1,4) = AsteroidH;
asteroidCloseApp(1,5) = 0;

end