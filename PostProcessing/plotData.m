function [] = plotData(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff,fuelMass)

%% Data that is required

au=1.49597870691E11; % conversion of m/au

% Solving differential motions
timeFinal=(6.653820100923719e+07); % orbital period
tspan=[timeFinal 0];
options = odeset('RelTol',1e-12);
[tE, yE] = ode45(@orbitalMotion,tspan,y0E,options,gammaCoeff,tauCoeff,timeFinal,0);
[tA, yA] = ode45(@orbitalMotion,tspan,y0A,options,gammaCoeff,tauCoeff,timeFinal,0);

% Transform to cartesian coordinates for position and velocity of asteroid and earth
[cX,cY,cZ]= pol2cart(cR(2,1:sizeC),cR(1,1:sizeC),cR(3,1:sizeC));
[eX,eY,eZ]= pol2cart(yE(:,2),yE(:,1),yE(:,3));
[aX,aY,aZ]= pol2cart(yA(:,2),yA(:,1),yA(:,3));
% Acceleration vector in cartesian coordinates
[accelX,accelY,accelZ] = getAccel(cR,tripTime,gammaCoeff,tauCoeff,sizeC);



%% Sub Plot 1

figure(1) %orbitals
subplot(2,3,1)
polarplot(yE(:,2),yE(:,1))
hold on
polarplot(yA(:,2),yA(:,1),'.')
hold on
polarplot(cR(2,1),cR(1,1),'r*')
hold on
polarplot(y0A(2),y0A(1),'*b')
hold on
polarplot(cR(2,1:sizeC),cR(1,1:sizeC))
legend({'earth','asteroid','launch','impact','spacecraft'})
title('r-\theta plane')
hold off

%radius vs. time
subplot(2,3,4)
plot(tE-(timeFinal-tripTime),yE(:,1))
hold on
plot(tA-(timeFinal-tripTime),yA(:,1),'.')
hold on
plot(cR(7,1:sizeC),cR(1,1:sizeC))
ylabel('r (a.u.)')
xlabel('t (s)')
xlim([0 tripTime])
title('Orbital radius')
hold off

%specific angular momentum vs. time
subplot(2,3,2)
plot(tE-(timeFinal-tripTime),yE(:,1).*yE(:,5))
hold on
plot(tA-(timeFinal-tripTime),yA(:,1).*yA(:,5),'.')
hold on
plot(cR(7,1:sizeC),cR(1,1:sizeC).*cR(5,1:sizeC))
ylabel('h (a.u.^{2}/s)')
xlabel('t (s)')
xlim([0 tripTime])
title('Specific angular momentum')
hold off

%plot z
subplot(2,3,3)
plot(yE(:,1).*cos(yE(:,2)),yE(:,3))
hold on
plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'.')
hold on
plot(cR(1,1:sizeC).*cos(cR(2,1:sizeC)), cR(3,1:sizeC),'LineWidth', 2)
xlabel('x (a.u.)')
ylabel('z (a.u.)')
title('x-z plane')
hold off

%Z vs. time
subplot(2,3,6)
plot(tE-(timeFinal-tripTime),yE(:,3))
hold on
plot(tA-(timeFinal-tripTime),yA(:,3),'.')
hold on
plot(cR(7,1:sizeC),cR(3,1:sizeC))
ylabel('z (a.u.)')
xlim([0 tripTime])
xlabel('t (s)')
title('Orbital elevation')
hold off

%theta vs. time
subplot(2,3,5)
plot(tE-(timeFinal-tripTime),mod(yE(:,2),2*pi))
hold on
plot(tA-(timeFinal-tripTime),mod(yA(:,2), 2*pi),'.')
hold on
plot(cR(7,1:sizeC),mod(cR(2,1:sizeC), 2*pi))
ylabel('\theta (rad.)')
xlim([0 tripTime])
xlabel('t (s)')
title('Orbital angle')
hold off

%% Subplot 2

figure(2) % acceleration vs. time
subplot(2,3,1)
plot(cR(7,1:sizeC),au*cR(10,1:sizeC))
xlim([0 tripTime])
title('Acceleration due to thrust')
ylabel('a_{thrust} (m/s^{2})')
xlabel('t (s)')

hold off

subplot(2,3,4)
plot(cR(7,1:sizeC),sin(cR(8,1:sizeC)).*cos(cR(9,1:sizeC)))
xlim([0 tripTime]), ylim([-1,1])
%legend('matlab','c')
title('Radial thrust fraction')
xlabel('t (s)')
ylabel('sin(\gamma)cos(\tau)')

subplot(2,3,5)
plot(cR(7,1:sizeC),cos(cR(8,1:sizeC)).*cos(cR(9,1:sizeC)))
xlim([0 tripTime]), ylim([-1,1])
%legend('matlab','c')
title('Tangential thrust fraction')
xlabel('t (s)')
ylabel('cos(\gamma)cos(\tau)')

subplot(2,3,6)
plot(cR(7,1:sizeC),sin(cR(9,1:sizeC)))
xlim([0 tripTime]), ylim([-1,1])
%legend('matlab','c')
title('Off-plane thrust fraction')
xlabel('t (s)')
ylabel('sin(\tau)')

%% coasting plots

% figure
% %plot(cR(7,1:sizeC),cR(11,1:sizeC))
% xlabel('Time')
% ylabel('coast')

% figure
% %5plot(cR(7,1:sizeC),cR(12,1:sizeC))
% xlabel('Time')
% ylabel('coast value')

co = angles(cR(7,1:sizeC),tripTime,coast);
subplot(2,3,3)
plot(cR(7,1:sizeC),sin(co).^2)
xlim([0 tripTime]), ylim([0,1])
title('Coasting function and threshold')
xlabel('t (s)')
ylabel('sin^2(\psi)')
hold on
coast_thresholdPlot = coast_threshold*ones(1,sizeC); % creates a vector with values of coast_threshold so MATLAB can plot it as a line
plot(cR(7,1:sizeC),coast_thresholdPlot,'--','color','r')
xlim([0 tripTime])
hold off

fuelSpent = (fuelMass - cR(11,1:sizeC))/fuelMass;
subplot(2,3,2)
plot(cR(7,1:sizeC),fuelSpent*100)
xlim([0 tripTime])
ylim([0,100])
title('Fuel consumption')
ylabel('% fuel (kg)')
xlabel('t (s)')

% Thrust angle plots
figure(3)

% Test Fourier calculation for start
% vpa(cR(8,1)), vpa(cR(9,1))

subplot(3,1,1)
plot(cR(7,1:sizeC),mod(cR(8,1:sizeC),2*pi))
xlabel('t (s)'), ylabel('\gamma (rad.)')
xlim([0 tripTime])
title('In-plane thrust angle')

subplot(3,1,2)
plot(cR(7,1:sizeC),cR(9,1:sizeC))
xlabel('t (s)'), ylabel('\tau (rad.)')
xlim([0 tripTime])
title('Out-of-plane thrust angle')

subplot(3,1,3)
plot(cR(7,1:sizeC),co)
xlabel('t (s)'), ylabel('\psi')
xlim([0 tripTime])
title('Coast series')

%% full orbital plots (vectors and no vectors)

radStep=1:15:length(cX)*1.0;
%a=figure(3); % plot with vectors
figure(4) % plot with vectors
plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880]	)
%plot3(cX(1),cY(1),cZ(1),'*','LineWidth', 5,'Color',[0.9290, 0.6940, 0.1250])
%hold on
%plot3(cX(end),cY(end),cZ(end),'*','LineWidth', 5,'Color',[0.9290, 0.6940, 0.1250])
xlim([-2.5 1.5])
ylim([-2.0 2.0])
zlim([-0.2 0.2])
xlabel('x (a.u.)')
ylabel('y (a.u.)')
zlabel('z (a.u.)')

hold on
plot3(aX,aY,aZ,'LineWidth', 1, 'Color',	[0.6350, 0.0780, 0.1840])
hold on
plot3(eX,eY,eZ,'LineWidth', 1,'Color',[.61 .51 .74])
hold on
quiver3(cX(radStep),cY(radStep),cZ(radStep),accelX(radStep),accelY(radStep),accelZ(radStep),'k','Autoscalefactor',.08,'LineWidth',1)
title('Solar orbitals')
legend('Spacecraft','Asteroid','Earth','Acceleration')
hold off
%print(a,'3D.png','-dpng','-r300'); 


%b=figure(4);
%figure(4)
%plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'LineWidth', 1,'Color',[.61 .51 .74]);
%xlim([-2.5 2])
%ylim([-.3 .3])
%hold on
%plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'LineWidth', 1,'Color',	[0.6350, 0.0780, 0.1840])
%hold on
%plot(cR(1,1:sizeC).*cos(cR(2,1:sizeC)), cR(3,1:sizeC),'LineWidth', 1,'Color',[0.4660, 0.6740, 0.1880])
%xlabel('x')
%ylabel('z')
%legend('Earth','Asteroid','Spacecraft')
%hold on
%quiver(cX(radStep),cZ(radStep),accelX(radStep),accelZ(radStep),'k','LineWidth',1,'MarkerSize',15,'Autoscalefactor',.08)
%title('Solar orbitals')
%hold off
%print(b,'2DNoVec.png','-dpng','-r350'); 

end

