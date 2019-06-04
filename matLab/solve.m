% Constants
Tperiod = 770.1180672365416; %Time period in days
sec=24*3600;
tspan=[0 Tperiod*sec];
au=1.49597870691E11;%Conversion between meters and AU

%Gets the X and Y coordinates of the asteroid
[Xa,Ya,Za,~,~,~] = ParseFile('AsteroidPosition','asteroid_');
%Gets the X and Y velocities of the asteroid
[Vxa,Vya,Vza,~,~,~] = ParseFile('AsteroidVelocity','asteroid_');

%Gets the X and Y coordinates of the earth
[Xe,Ye,Ze,~,~,~] = ParseFile('EarthPosition','earth_');
%Gets the X and Y velocities of the earth
[Vxe,Vye,Vze,~,~,~] = ParseFile('EarthVelocity','earth_');

% Position transformation for asteroid
xCartesian=[Xa(1);Ya(1);Za(1)]*au; % Cartesian positions for asteroid
[theta1, rho1, z1] = cart2pol(xCartesian(1), xCartesian(2), xCartesian(3));
% Velocity transformation for asteroid
vCartesian=[Vxa(1);Vya(1);Vza(1)]*au/sec; % Cartesian velocities for asteroid
% Transformation to cylindrical for asteroid
vMatrix=[cos(theta1)  sin(theta1) 0
         -sin(theta1) cos(theta1) 0
         0           0          1];
% Velocities in cylindrical coordinates for asteroid
V=vMatrix*vCartesian;

% Position transformation for earth
xCartesianE=[Xe(1);Ye(1);Ze(1)]*au; % Cartesian positions for earth
[thetaE, rhoE, zE] = cart2pol(xCartesianE(1), xCartesianE(2), xCartesianE(3));
% Velocity transformation for earth
vCartesianE=[Vxe(1);Vye(1);Vze(1)]*au/sec; % Cartesian velocities for earth
% Transformation to cylindrical for earth
vMatrixE=[cos(thetaE)  sin(thetaE) 0
         -sin(thetaE) cos(thetaE) 0
         0           0          1];
% Velocities in cylindrical coordinates for earth
Ve=vMatrixE*vCartesianE;
     
% asteroid plot data
% [theta,rho] = cart2pol(Xa*au,Ya*au);
% [vtheta,vrho] = cart2pol(Vxa*au/sec,Vya*au/sec);
% Initial conditions for asteroid
% y0=[rho1 theta1 z1 V(1) V(2) V(3)];

% earth plot data
[theta,rho] = cart2pol(Xe*au,Ye*au);
[vtheta,vrho] = cart2pol(Vxe*au/sec,Vye*au/sec);
% Initial conditions for earth
 y0=[rhoE(1) thetaE(1) Ze(1) Ve(1) Ve(2) Ve(3)];

%Solving differential motions
options = odeset('RelTol',1e-9);
tic
[t, y] = ode45(@orbitalMotion,tspan,y0, options);
toc

%figures
figure(1) %orbitals
polarplot(y(:,2),y(:,1),'--')
hold on
polarplot(theta,rho)
hold off

% figure(3) % Conservation of momentum
% h=y(:,1).*y(:,5);
% plot(t,h)