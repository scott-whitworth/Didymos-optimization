function [] = plotDataCompare(cR,y0A,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,fuelMass1,alpha1,beta1,zeta1,launchPos1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2,fuelMass2,alpha2,beta2,zeta2,launchPos2)
    %% Data that is required
    
    au=1.49597870691E11; % conversion of m/au
    
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
    polarplot(yE(:,2),yE(:,1),'.')
    hold on
    polarplot(yA(:,2),yA(:,1),'.')
    hold on
    polarplot(cR(2,1),cR(1,1),'r*')
    hold on
    polarplot(dR(2,1),dR(1,1),'g*')
    hold on
    polarplot(y0A(2),y0A(1),'*b')
    hold on
    polarplot(cR(2,:),cR(1,:))
    hold on
    polarplot(dR(2,:),dR(1,:))
    legend({'earth','asteroid','launch 1','launch 2','impact','spacecraft 1','spacecraft 2'})
    title('r-\theta plane')
    hold off
    
    %radius vs. time
    minTripTime = min(tripTime1,tripTime2);
    maxTripTime = max(tripTime1,tripTime2);
    tripTimeDiff = maxTripTime - minTripTime;
    subplot(2,3,4)
    plot(tE-(timeFinal-maxTripTime),yE(:,1),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,1),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(1,:))
    else
        plot(cR(7,:),cR(1,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(1,:))
    else
        plot(dR(7,:),dR(1,:))
    end
    ylabel('r (a.u.)')
    xlabel('t (s)')
    xlim([0 maxTripTime])
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital radius')
    hold off
    
    %angular momentum vs. time
    subplot(2,3,2)
    plot(tE-(timeFinal-maxTripTime),yE(:,1).*yE(:,5),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,1).*yA(:,5),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(1,:).*cR(5,:))
    else
        plot(cR(7,:),cR(1,:).*cR(5,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(1,:).*dR(5,:))
    else
        plot(dR(7,:),dR(1,:).*dR(5,:))
    end
    ylabel('h (a.u.^{2}/s)')
    xlabel('t (s)')
    xlim([0 maxTripTime])
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Specific Angular Momentum')
    hold off
    
    %plot z
    subplot(2,3,3)
    plot(yE(:,1).*cos(yE(:,2)),yE(:,3),'.')
    hold on
    plot(yA(:,1).*cos(yA(:,2)),yA(:,3),'.')
    hold on
    plot(cR(1,:).*cos(cR(2,:)), cR(3,:),'LineWidth', 2)
    hold on
    plot(dR(1,:).*cos(dR(2,:)), dR(3,:),'LineWidth', 2)
    xlabel('x (a.u.)')
    ylabel('z (a.u.)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('x-z plane')
    hold off
    
    %Z vs. time
    subplot(2,3,6)
    plot(tE-(timeFinal-maxTripTime),yE(:,3),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),yA(:,3),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,cR(3,:))
    else
        plot(cR(7,:),cR(3,:))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,dR(3,:))
    else
        plot(dR(7,:),dR(3,:))
    end
    ylabel('z (a.u.)')
    xlim([0 maxTripTime])
    xlabel('t (s)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital elevation')
    hold off
    
    %theta vs. time
    subplot(2,3,5)
    plot(tE-(timeFinal-maxTripTime),mod(yE(:,2),2*pi),'.')
    hold on
    plot(tA-(timeFinal-maxTripTime),mod(yA(:,2), 2*pi),'.')
    hold on
    if tripTime1 < tripTime2
        plot(cR(7,:)+tripTimeDiff,mod(cR(2,:), 2*pi))
    else
        plot(cR(7,:),mod(cR(2,:), 2*pi))
    end
    hold on
    if tripTime2 < tripTime1
        plot(dR(7,:)+tripTimeDiff,mod(dR(2,:), 2*pi))
    else
        plot(dR(7,:),mod(dR(2,:), 2*pi))
    end
    ylabel('\theta (rad.)')
    xlim([0 maxTripTime])
    xlabel('t (s)')
    legend({'earth','asteroid','spacecraft 1','spacecraft 2'})
    title('Orbital angle')
    hold off
    
    %% Subplot 2
    
    figure(2) % only spacecraft 1
    subplot(2,2,1)
    plot(cR(7,:),au*cR(10,:))
    xlim([0 tripTime1])
    legend({'spacecraft 1'})
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2}')
    xlabel('t (s)')
    
    hold off
    
    
    co = angles(cR(7,:),tripTime1,coast1);
    subplot(2,2,2)
    plot(cR(7,:),sin(co).^2)
    title('Coasting function and threshold')
    xlabel('t (s)')
    ylabel('sin^2(\psi)')
    hold on
    coast_threshold1Plot = coast_threshold1*ones(1,sizeC);
    plot(cR(7,:),coast_threshold1Plot,'--','color','r')
    xlim([0 tripTime1]), ylim([0,1])
    legend({'spacecraft 1','threshold 1'})
    hold off
    
    fuelSpent1 = (fuelMass1 - cR(11,:))/fuelMass1;
    subplot(2,2,3)
    plot(cR(7,:),fuelSpent1*100)
    xlim([0 tripTime1])
    title('Fuel consumption')
    legend({'spacecraft 1'})
    ylabel('% fuel')
    xlabel('t (s)')

    err1 = (cR(12,:)-cR(13,:))./cR(14,:);
    subplot(2,2,4)
    plot(cR(7,:),err1*100)
    xlim([0 tripTime1])
    title('Conservation of mechanical energy')
    legend({'spacecraft 1'})
    ylabel('% error')
    xlabel('t (s)')
    
    figure(3) % only spacecraft 2
    subplot(2,2,1)
    plot(dR(7,:),au*dR(10,:))
    xlim([0 tripTime2])
    legend({'spacecraft 2'})
    title('Acceleration due to thrust')
    ylabel('a_{thrust} (m/s^{2})')
    xlabel('t (s)')
    
    hold off
    
    
    do = angles(dR(7,:),tripTime2,coast2);
    subplot(2,2,2)
    plot(dR(7,:),sin(do).^2)
    title('Coasting function and threshold')
    xlabel('t (s)')
    ylabel('sin^2(\psi)')
    hold on
    coast_threshold2Plot = coast_threshold2*ones(1,sizeD);
    plot(dR(7,:),coast_threshold2Plot,'--','color','r')
    xlim([0 tripTime2]), ylim([0,1])
    legend({'spacecraft 2','threshold 2'})
    hold off
    
    fuelSpent2 = (fuelMass2 - dR(11,:))/fuelMass2;
    subplot(2,2,3)
    plot(dR(7,:),fuelSpent2*100)
    xlim([0 tripTime2])
    title('Fuel consumption')
    legend({'spacecraft 2'})
    ylabel('% fuel')
    xlabel('t (s)')

    err2 = (dR(12,:)-dR(13,:))./dR(14,:);
    subplot(2,2,4)
    plot(dR(7,:),err2*100)
    xlim([0 tripTime2])
    title('Conservation of mechanical energy')
    legend({'spacecraft 2'})
    ylabel('% error')
    xlabel('t (s)')

    % only spacecraft 1
    figure(4)

    subplot(2,3,1)
    plot(cR(7,:), au*cR(4,:))
    xlim([0 tripTime1])
    title('Radial velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{r} (m/s)')

    subplot(2,3,2)
    plot(cR(7,:), au*cR(5,:))
    xlim([0 tripTime1])
    title('Tangential velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{\theta} (m/s)')

    subplot(2,3,3)
    plot(cR(7,:), au*cR(6,:))
    xlim([0 tripTime1])
    title('Axial velocity')
    legend({'spacecraft 1'})
    xlabel('t (s)')
    ylabel('v_{z} (m/s)')
    
    subplot(2,3,4)
    plot(cR(7,:),sin(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Radial thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\gamma)cos(\tau)')
    
    subplot(2,3,5)
    plot(cR(7,:),cos(cR(8,:)).*cos(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Tangential thrust fraction')
    xlabel('t (s)')
    ylabel('cos(\gamma)cos(\tau)')
    
    subplot(2,3,6)
    plot(cR(7,:),sin(cR(9,:)))
    xlim([0 tripTime1]), ylim([-1,1])
    legend({'spacecraft 1'})
    title('Off-plane thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\tau)')
    
    % only spacecraft 2
    figure(5)

    subplot(2,3,1)
    plot(dR(7,:), au*dR(4,:))
    xlim([0 tripTime2])
    title('Radial velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{r} (m/s)')

    subplot(2,3,2)
    plot(dR(7,:), au*dR(5,:))
    xlim([0 tripTime2])
    title('Tangential velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{\theta} (m/s)')

    subplot(2,3,3)
    plot(dR(7,:), au*dR(6,:))
    xlim([0 tripTime2])
    title('Axial velocity')
    legend({'spacecraft 2'})
    xlabel('t (s)')
    ylabel('v_{z} (m/s)')
    
    subplot(2,3,4)
    plot(dR(7,:),sin(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Radial thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\gamma)cos(\tau)')
    
    subplot(2,3,5)
    plot(dR(7,:),cos(dR(8,:)).*cos(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Tangential thrust fraction')
    xlabel('t (s)')
    ylabel('cos(\gamma)cos(\tau)')
    
    subplot(2,3,6)
    plot(dR(7,:),sin(dR(9,:)))
    xlim([0 tripTime2]), ylim([-1,1])
    legend({'spacecraft 2'})
    title('Off-plane thrust fraction')
    xlabel('t (s)')
    ylabel('sin(\tau)')

    % Thrust angle plots
    figure(6)
    
    % Gamma plots
    subplot(3,1,1)
    plot(cR(7,:),mod(cR(8,:),2*pi))
    hold on
    plot(dR(7,:),mod(dR(8,:),2*pi))
    xlabel('t (s)'), ylabel('\gamma (rad.)')
    xlim([0 maxTripTime])
    legend('spacecraft 1', 'spacecraft 2')
    title('In-plane thrust angle')
    hold off
    
    % Tau plots
    subplot(3,1,2)
    plot(cR(7,:),cR(9,:))
    hold on
    plot(dR(7,:),dR(9,:))
    xlabel('t (s)'), ylabel('\tau (rad.)')
    xlim([0 maxTripTime])
    legend('spacecraft 1', 'spacecraft 2')
    title('Out-of-plane thrust angle')
    hold off
    
    % Psi plots
    subplot(3,1,3)
    plot(cR(7,:),co)
    hold on
    plot(dR(7,:),do)
    xlabel('t (s)'), ylabel('\psi')
    xlim([0 maxTripTime])
    legend('spacecraft 1', 'spacecraft 2')
    title('Coast series')
    hold off
    
    %% full orbital plots (vectors and no vectors)
    
    radStep1=1:15:length(cX)*1.0;
    radStep2=1:15:length(dX)*1.0;
    %a=figure(3); % plot with vectors
    figure(7) % plot with vectors
    plot3(cX,cY,cZ,'LineWidth', 3,'Color',[0.4660, 0.6740, 0.1880]	)
    hold on
    plot3(dX,dY,dZ,'LineWidth', 3,'Color',[0, 0, 1]	)
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

    figure(8)
    
    r_esoi = 6.211174738e-3; % radius of Earth's sphere of influence in au
    t = 0:pi/100:2*pi;
    [m,n] = size(t);
    
    plot3(launchPos1(1)+r_esoi*cos(t), launchPos1(2)+r_esoi*sin(t), launchPos1(3)*ones(1,n),'LineWidth',1,'Color',[1 0 0])
    hold on
    [alpha_x1, alpha_y1, alpha_z1] = pol2cart(alpha1, r_esoi, 0);
    plot3(alpha_x1+launchPos1(1), alpha_y1+launchPos1(2), alpha_z1+launchPos1(3),'*','Color',[1 0 0])
    hold on
    quiver3(alpha_x1+launchPos1(1), alpha_y1+launchPos1(2), alpha_z1+launchPos1(3), sin(beta1)*cos(zeta1), cos(beta1)*cos(zeta1), sin(zeta1),'Autoscalefactor',.005,'LineWidth',1,'Color',[1 0 0]);
    hold on
    
    plot3(launchPos2(1)+r_esoi*cos(t), launchPos2(2)+r_esoi*sin(t), launchPos2(3)*ones(1,n),'LineWidth',1,'Color',[0 0 1])
    hold on
    [alpha_x2, alpha_y2, alpha_z2] = pol2cart(alpha2, r_esoi, 0);
    plot3(alpha_x2+launchPos2(1), alpha_y2+launchPos2(2), alpha_z2+launchPos2(3),'*','Color',[0 0 1])
    hold on
    quiver3(alpha_x2+launchPos2(1), alpha_y2+launchPos2(2), alpha_z2+launchPos2(3), sin(beta2)*cos(zeta2), cos(beta2)*cos(zeta2), sin(zeta2),'Autoscalefactor',.005,'LineWidth',1,'Color',[0 0 1]);
    hold on
   
    max_x = launchPos1(1); max_y = launchPos1(2); max_z = launchPos1(3);
    min_x = launchPos2(1); min_y = launchPos2(2); min_z = launchPos2(3);
    if max_x < min_x
        max_x = launchPos2(1);
        min_x = launchPos1(1);
    end
    if max_y < min_y
        max_y = launchPos2(2);
        min_y = launchPos1(2);
    end
    if max_z < min_z
        max_z = launchPos2(3);
        min_z = launchPos1(3);
    end
    diff_x = max_x - min_x; diff_y = max_y - min_y; diff_z = max_z - min_z;
    if diff_x > diff_y && diff_x > diff_z
        diff_max = diff_x;
    elseif diff_y > diff_z
        diff_max = diff_y;
    else
        diff_max = diff_z;
    end
    scale_x = (diff_max-diff_x)/2 + 2*r_esoi; scale_y = (diff_max-diff_y)/2 + 2*r_esoi; scale_z = (diff_max-diff_z)/2 + 2*r_esoi;
    
    % xlim([min_x-2*r_esoi, max_x+2*r_esoi])
    % ylim([min_y-2*r_esoi, max_y+2*r_esoi])
    % zlim([min_z-2*r_esoi, max_z+2*r_esoi]) 
    xlim([min_x-scale_x, max_x+scale_x])
    ylim([min_y-scale_y, max_y+scale_y])
    zlim([min_z-scale_z, max_z+scale_z])
    xlabel('x (a.u.)')
    ylabel('y (a.u.)')
    zlabel('z (a.u.)')
    title('Launch conditions')
    legend({'ESOI_1','Position_1','Velocity_1','ESOI_2','Position_2','Velocity_2'})
    hold on
    
    end