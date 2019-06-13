clear ALL
%% Load file

fileID = fopen('orbitalMotion-accel.bin');
sizeC=966;
cR = fread(fileID,[9 sizeC],'double');

%% Constants
Tperiod = 770.1180672365416; %Time period in days
sec=24*3600;
timeFinal=6.653820100923719e+07;
tspan=[0 timeFinal];
au=1.49597870691E11;%Conversion between meters and AU

%% Initial conditions
%Gets the X and Y coordinates of the asteroid
[Xa,Ya,Za,~,~,~] = ParseFile('AsteroidPosition','asteroid_');
%Gets the X and Y velocities of the asteroid
[Vxa,Vya,Vza,~,~,~] = ParseFile('AsteroidVelocity','asteroid_');

%Gets the X and Y coordinates of the earth
[Xe,Ye,Ze,~,~,~] = ParseFile('EarthPosition','earth_');
%Gets the X and Y velocities of the earth
[Vxe,Vye,Vze,~,~,~] = ParseFile('EarthVelocity','earth_');

[R,V] = polarT(Xe(1),Ye(1),Ze(1),Vxe(1),Vye(1),Vze(1));
y0E = [R V/sec];
[R,V] = polarT(Xa(1),Ya(1),Za(1),Vxa(1),Vya(1),Vza(1));
%y0A = [R V/sec];
y0A = [1.02696822710421, 0.238839574416454, -0.0526614832914496,...
    -2.05295246185041e-08, 2.29132593453064e-07, 8.00663905822009e-09];
y0S = [cR(1,1), cR(2,1), cR(3,1),...
   cR(4,1), cR(5,1), cR(6,1)];

V_ini = sqrt(y0S(4)^2+y0S(5)^2+y0S(6)^2);
V_ini = V_ini*au;
disp(V_ini)

%% Solving differential motions

gammaCoeff= [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
tauCoeff=[1,1,1,1,1,1,1,1,1,1,1];
accel=0.000/au;
options = odeset('RelTol',1e-15);
tic
[tS, yS] = ode45(@orbitalMotion,tspan,y0S, options,gammaCoeff,tauCoeff,timeFinal,accel);
[tE, yE] = ode45(@orbitalMotion,tspan,y0E, options,gammaCoeff,tauCoeff,timeFinal,0);
[tA, yA] = ode45(@orbitalMotion,tspan,y0A, options,gammaCoeff,tauCoeff,timeFinal,0);
toc

%% Figures
figure(1) %orbitals
polarplot(yS(:,2),yS(:,1),'--')
hold on
polarplot(yE(:,2),yE(:,1))
hold on
polarplot(yA(:,2),yA(:,1),'.')
hold on
polarplot(.238839574416454,1.02696822710421,'*r')
hold on
polarplot(cR(2,1:sizeC),cR(1,1:sizeC))
legend({'spaceCraft','earth','asteroid'})
hold off

figure(2)
plot(tS,yS(:,1),'--')
hold on
plot(tE,yE(:,1))
hold on
plot(tA,yA(:,1),'.')
hold on
plot(cR(7,1:sizeC),cR(1,1:sizeC))
hold on
plot(tS(end),1.02696822710421,'r*')
legend({'spaceCraft','earth','asteroid','sC-c'})
ylabel('r')
xlabel('t')
hold off

figure(3)
plot(tS,yS(:,1).*yS(:,5),'--')
hold on
plot(tE,yE(:,1).*yE(:,5))
hold on
plot(tA,yA(:,1).*yA(:,5),'.')
hold on
plot(cR(7,1:sizeC),cR(1,1:sizeC).*cR(5,1:sizeC))
legend({'spaceCraft','earth','asteroid','sC-c'})
ylabel('h')
xlabel('t')
hold off


%% Angles

% [ga,ta]=angles(tS,Tperiod*sec/2,gammaCoeff,tauCoeff);
% figure (4)
% plot(tS,ga)
% legend('matlab','c')
% title('gamma')
% xlabel('Time')
% ylabel('gamma')
% 
% figure (5)
% plot(tS,ta)
% legend('matlab','c')
% title('tau')
% xlabel('Time')
% ylabel('tau')

