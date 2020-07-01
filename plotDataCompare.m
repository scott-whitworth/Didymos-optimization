function [] = plotDataCompare(cR,y0A,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2)
%% Data that is required

% Solving differential motions
timeFinal=(6.653820100923719e+07);
tspan=[timeFinal 0];
options = odeset('RelTol',1e-12);
[tE, yE] = ode45(@orbitalMotion,tspan,y0E,options,gammaCoeff1,tauCoeff1,timeFinal,0);
[tA, yA] = ode45(@orbitalMotion,tspan,y0A,options,gammaCoeff1,tauCoeff1,timeFinal,0);

% Transform to cartesian coordinates for position and velocity of asteroid and earth
[cX,cY,cZ]= pol2cart(cR(2,1:sizeC),cR(1,1:sizeC),cR(3,1:sizeC));
[dX,dY,dZ]= pol2cart(dR(2,1:sizeD),dR(1,1:sizeD),dR(3,1:sizeD));
[eX,eY,eZ]= pol2cart(yE(:,2),yE(:,1),yE(:,3));
[aX,aY,aZ]= pol2cart(yA(:,2),yA(:,1),yA(:,3));
% Acceleration vector in cartesian coordinates
[accelX1,accelY1,accelZ1] = getAccel(cR,tripTime1,gammaCoeff1,tauCoeff1,sizeC);
[accelX2,accelY2,accelZ2] = getAccel(dR,tripTime2,gammaCoeff2,tauCoeff2,sizeD);



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
hold on
polarplot(dR(2,1:sizeD),dR(1,1:sizeD))
legend({'earth','asteroid','launch','impact','spacecraft 1','spacecraft 2'})
title('r-\theta plane')
hold off

%radius vs. time
minTripTime = min(tripTime1,tripTime2);
maxTripTime = max(tripTime1,tripTime2);
tripTimeDiff = maxTripTime - minTripTime;
subplot(2,3,4)
plot(tE-(timeFinal-maxTripTime),yE(:,1))
hold on
plot(tA-(timeFinal-maxTripTime),yA(:,1),'.')
hold on
if tripTime1 < tripTime2
    plot(cR(7,1:sizeC)+tripTimeDiff,cR(1,1:sizeC))
else
    plot(cR(7,1:sizeC),cR(1,1:sizeC))
end
hold on
if tripTime2 < tripTime1
    plot(dR(7,1:sizeD)+tripTimeDiff,dR(1,1:sizeD))
else
    plot(dR(7,1:sizeD),dR(1,1:sizeD))
end
ylabel('r')
xlabel('t')
xlim([0 maxTripTime])
legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
title('Orbital radius')
hold off

%angular momentum vs. time
subplot(2,3,2)
plot(tE-(timeFinal-maxTripTime),yE(:,1).*yE(:,5))
hold on
plot(tA-(timeFinal-maxTripTime),yA(:,1).*yA(:,5),'.')
hold on
if tripTime1 < tripTime2
    plot(cR(7,1:sizeC)+tripTimeDiff,cR(1,1:sizeC).*cR(5,1:sizeC))
else
    plot(cR(7,1:sizeC),cR(1,1:sizeC).*cR(5,1:sizeC))
end
hold on
if tripTime2 < tripTime1
    plot(dR(7,1:sizeD)+tripTimeDiff,dR(1,1:sizeD).*dR(5,1:sizeD))
else
    plot(dR(7,1:sizeD),dR(1,1:sizeD).*dR(5,1:sizeD))
end
ylabel('L')
xlabel('t')
xlim([0 maxTripTime])
legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
title('Specific Angular Momentum')
hold off

%plot z
subplot(2,3,3)
plot(yE(:,1).*cos(yE(:,2)),yE(:,3))
hold on
plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'.')
hold on
plot(cR(1,1:sizeC).*cos(cR(2,1:sizeC)), cR(3,1:sizeC),'LineWidth', 2)
hold on
plot(dR(1,1:sizeD).*cos(dR(2,1:sizeD)), dR(3,1:sizeD),'LineWidth', 2)
xlabel('x')
ylabel('z')
legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
title('x-z plane')
hold off

%Z vs. time
subplot(2,3,6)
plot(tE-(timeFinal-maxTripTime),yE(:,3))
hold on
plot(tA-(timeFinal-maxTripTime),yA(:,3),'.')
hold on
if tripTime1 < tripTime2
    plot(cR(7,1:sizeC)+tripTimeDiff,cR(3,1:sizeC))
else
    plot(cR(7,1:sizeC),cR(3,1:sizeC))
end
hold on
if tripTime2 < tripTime1
    plot(dR(7,1:sizeD)+tripTimeDiff,dR(3,1:sizeD))
else
    plot(dR(7,1:sizeD),dR(3,1:sizeD))
end
ylabel('z')
xlim([0 maxTripTime])
xlabel('t')
legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
title('Orbital elevation')
hold off

%theta vs. time
subplot(2,3,5)
plot(tE-(timeFinal-maxTripTime),mod(yE(:,2),2*pi))
hold on
plot(tA-(timeFinal-maxTripTime),mod(yA(:,2), 2*pi),'.')
hold on
if tripTime1 < tripTime2
    plot(cR(7,1:sizeC)+tripTimeDiff,mod(cR(2,1:sizeC), 2*pi))
else
    plot(cR(7,1:sizeC),mod(cR(2,1:sizeC), 2*pi))
end
hold on
if tripTime2 < tripTime1
    plot(dR(7,1:sizeD)+tripTimeDiff,mod(dR(2,1:sizeD), 2*pi))
else
    plot(dR(7,1:sizeD),mod(dR(2,1:sizeD), 2*pi))
end
ylabel('\theta')
xlim([0 maxTripTime])
xlabel('t')
legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
title('Orbital angle')
hold off

%% Subplot 2

figure(2) % only spacecraft 1
subplot(2,3,1)
plot(cR(7,1:sizeC),cR(10,1:sizeC))
xlim([0 tripTime1])
legend({'spacecraft 1'})
title('Acceleration due to thrust')
ylabel('a_{thrust}')
xlabel('t')
hold off

subplot(2,3,4)
plot(cR(7,1:sizeC),sin(cR(8,1:sizeC)).*cos(cR(9,1:sizeC)))
xlim([0 tripTime1])
legend({'spacecraft 1'})
title('Radial thrust fraction')
xlabel('t')
ylabel('sin(\gamma)cos(\tau)')

subplot(2,3,5)
plot(cR(7,1:sizeC),cos(cR(8,1:sizeC)).*cos(cR(9,1:sizeC)))
xlim([0 tripTime1])
legend({'spacecraft 1'})
title('Tangential thrust fraction')
xlabel('t')
ylabel('cos(\gamma)cos(\tau)')

subplot(2,3,6)
plot(cR(7,1:sizeC),sin(cR(9,1:sizeC)))
xlim([0 tripTime1])
legend({'spacecraft 1'})
title('Off-plane thrust fraction')
xlabel('t')
ylabel('sin(\tau)')


[co,~]=angles(cR(7,1:sizeC),tripTime1,coast1,coast1);
subplot(2,3,2:3)
plot(cR(7,1:sizeC),sin(co).^2)
title('Coasting function and threshold')
xlabel('t')
ylabel('sin^2(\psi)')
hold on
plot(cR(7,1:sizeC),coast_threshold1)
xlim([0 tripTime1])
legend({'spacecraft 1','threshold 1'})
hold off

figure(3) % only spacecraft 2
subplot(2,3,1)
plot(dR(7,1:sizeD),dR(10,1:sizeD))
xlim([0 tripTime2])
legend({'spacecraft 2'})
title('Acceleration due to thrust')
ylabel('a_{thrust}')
xlabel('t')
hold off

subplot(2,3,4)
plot(dR(7,1:sizeD),sin(dR(8,1:sizeD)).*cos(dR(9,1:sizeD)))
xlim([0 tripTime2])
legend({'spacecraft 2'})
title('Radial thrust fraction')
xlabel('t')
ylabel('sin(\gamma)cos(\tau)')

subplot(2,3,5)
plot(dR(7,1:sizeD),cos(dR(8,1:sizeD)).*cos(dR(9,1:sizeD)))
xlim([0 tripTime2])
legend({'spacecraft 2'})
title('Tangential thrust fraction')
xlabel('t')
ylabel('cos(\gamma)cos(\tau)')

subplot(2,3,6)
plot(dR(7,1:sizeD),sin(dR(9,1:sizeD)))
xlim([0 tripTime2])
legend({'spacecraft 2'})
title('Off-plane thrust fraction')
xlabel('t')
ylabel('sin(\tau)')


[do,~]=angles(dR(7,1:sizeD),tripTime2,coast2,coast2);
subplot(2,3,2:3)
plot(dR(7,1:sizeD),sin(do).^2)
title('Coasting function and threshold')
xlabel('t')
ylabel('sin^2(\psi)')
hold on
plot(dR(7,1:sizeD),coast_threshold2)
xlim([0 tripTime2])
legend({'spacecraft 2','threshold 2'})
hold off


%% full orbital plots (vectors and no vectors)

radStep1=1:15:length(cX)*1.0;
radStep2=1:15:length(dX)*1.0;
%a=figure(3); % plot with vectors
figure(4) % plot with vectors
plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880]	)
hold on
plot3(dX,dY,dZ,'LineWidth', 3,'Color',[0, 0, 1]	)
xlim([-2.5 1.5])
ylim([-2.0 2.0])
zlim([-0.2 0.2])
xlabel('x')
ylabel('y')
zlabel('z')

hold on
plot3(aX,aY,aZ,'LineWidth', 1, 'Color',	[0.6350, 0.0780, 0.1840])
hold on
plot3(eX,eY,eZ,'LineWidth', 1,'Color',[.61 .51 .74])
hold on
quiver3(cX(radStep1),cY(radStep1),cZ(radStep1),accelX1(radStep1),accelY1(radStep1),accelZ1(radStep1),'k','Autoscalefactor',.08,'LineWidth',1,'Color',[0.4660, 0.6740, 0.1880])
hold on
quiver3(dX(radStep2),dY(radStep2),dZ(radStep2),accelX2(radStep2),accelY2(radStep2),accelZ2(radStep2),'k','Autoscalefactor',.08,'LineWidth',1,'Color',[0, 0, 1])
title('Solar orbitals')
legend('Spacecraft 1','Spacecraft 2','Asteroid','Earth','Acceleration 1', 'Acceleration 2')
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

