%GETCLOSEAPPROACH Calculates index of closest approach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function GetCloseApproach takes table data for Earth and the asteroid's
%position and calculates the index at which closest approach occurs. The
%function then returns vectors containing the X and Y position and velocity
%data at this index for both the asteroid and Earth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by GetFinalState.
function [AsteroidValues,EarthValues,AsteroidEccentricity,AsteroidSemimajor] = GetCloseApproach(astDes)

%Gets the X and Y coordinates of the asteroid
[Xa,Ya,AsteroidEccentricity,AsteroidSemimajor,~] = ParseFile('AsteroidPosition',astDes);

%Gets the X and Y velocities of the asteroid
[Vxa,Vya,~,~,~] = ParseFile('AsteroidVelocity',astDes);

%Gets the X and Y coordinates of Earth
[Xe,Ye,~,~,~] = ParseFile('EarthPosition',astDes);

%Gets the X and Y velocities of Earth
[Vxe,Vye,~,~,~] = ParseFile('EarthVelocity',astDes);

%The usage of Xe is arbitrary, they are all the same size
for i = 1:1:size(Xe,2)
    %calculates distance between asteroid and Earth position at a point in
    %time with respect to the sun.
    s(i) = sqrt((Xa(i)-Xe(i))^2 + (Ya(i)-Ye(i))^2);
end

%Returns the index at which closest approach occurs
[~,index] = min(s);
    
%Converts X and Y at closest approach for Earth and asteroid to meters, then
%non-dimensionalizes them
N_Xa = Xa(index)*1E3/Constants.MCONVERSION;
N_Ya = Ya(index)*1E3/Constants.MCONVERSION;
N_Xe = Xe(index)*1E3/Constants.MCONVERSION;
N_Ye = Ye(index)*1E3/Constants.MCONVERSION;

%Converts X and Y velocity at closest approach for Earth and asteroid to meters per second, then
%non-dimensionalizes them
N_Vxa = Vxa(index)*1E3/Constants.MCONVERSION*Constants.SCONVERSION;
N_Vya = Vya(index)*1E3/Constants.MCONVERSION*Constants.SCONVERSION;
N_Vxe = Vxe(index)*1E3/Constants.MCONVERSION*Constants.SCONVERSION;
N_Vye = Vye(index)*1E3/Constants.MCONVERSION*Constants.SCONVERSION;

%Generates vectors containing X,Y,Vx,Vy for both the asteroid and Earth
AsteroidValues = [N_Xa,N_Ya,N_Vxa,N_Vya];
EarthValues = [N_Xe,N_Ye,N_Vxe,N_Vye];

end