function [] = individualProgress(str)
% individualProgress: plots generational data from a binary file
    % Input:
        % str: the specific generation to be plotted ('Best' or 'Worst')
    % Output:
        % N/A: see figures

fin = fopen(join([str,'InGenerations.bin']));
A = fread(fin, [20, Inf], 'double');

au = 1.49587870691e11;
km = au/1e3;

% Position graphs
figure(1)
% R initial over generations
sp1 = subplot(2,3,1);
plot(A(1,:), km*A(10,:))
title('Init. Radius')
xlabel('generations'), ylabel('r_i (km)')
% Theta initial over generations
sp2 = subplot(2,3,2);
plot(A(1,:), A(11,:))
title('Init. Direction')
xlabel('generations'), ylabel('\theta_i (rad.)')
% Z initial over generations
sp3 = subplot(2,3,3);
plot(A(1,:), km*A(12,:))
xlabel('generations'), ylabel('z_i (km)')
title('Init. Elevation')
% R final over generations
sp4 = subplot(2,3,4);
plot(A(1,:), km*A(4,:))
title('Final Radius')
xlabel('generations'), ylabel('r_f (km)')
% Theta final over generations
sp5 = subplot(2,3,5);
plot(A(1,:), A(5,:))
title('Final Direction')
xlabel('generations'), ylabel('\theta_f (rad.)')
% Z final over generations
sp6 = subplot(2,3,6);
plot(A(1,:), km*A(6,:))
title('Final Elevation')
xlabel('generations'), ylabel('z_f (km)')

linkaxes([sp1,sp2,sp3,sp4,sp5,sp6],'x')
suptitle('Initial (_i) and Final (_f) Positions')

% Velocity graphs
figure(2)
% VR initial over generations
vsp1 = subplot(2,3,1);
plot(A(1,:), km*A(13,:))
title('Init. Radial velocity')
xlabel('generations'), ylabel('v_r_i (km/s)')
% VTh initial over generations
vsp2 = subplot(2,3,2);
plot(A(1,:), km*A(14,:))
title('Init. Angular velocity')
xlabel('generations'), ylabel('v_\theta_i (km/s)')
% VZ initial over generations
vsp3 = subplot(2,3,3);
plot(A(1,:), km*A(15,:))
title('Init. Upward velocity')
xlabel('generations'), ylabel('v_z_i (km/s)')
% VR final over generations
vsp4 = subplot(2,3,4);
plot(A(1,:), km*A(7,:))
title('Final Radial velocity')
xlabel('generations'), ylabel('v_r_f (km/s)')
% VTh final over generations
vsp5 = subplot(2,3,5);
plot(A(1,:), km*A(8,:))
title('Final Angular velocity')
xlabel('generations'), ylabel('v_\theta_f (km/s)')
% VZ final over generations
vsp6 = subplot(2,3,6);
plot(A(1,:), km*A(9,:))
title('Final Upward velocity')
xlabel('generations'), ylabel('v_z_f (km/s)')

linkaxes([vsp1,vsp2,vsp3,vsp4,vsp5,vsp6],'x')
suptitle('Initial (_i) and Final (_f) Velocity')

% Initial parameter graphs
figure(3)
% Alpha over generations
psp1 = subplot(4,1,1);
plot(A(1,:), A(16,:))
title('Initial position angle')
xlabel('generations'), ylabel('\alpha_0 (rad.)')
% Beta over generations
psp2 = subplot(4,1,2);
plot(A(1,:), A(17,:))
title('Initial in-plane velocity angle')
xlabel('generations'), ylabel('\beta_0 (rad.)')
% Zeta over generations
psp3 = subplot(4,1,3);
plot(A(1,:), A(18,:))
title('Initial out-of-plane velocity angle')
xlabel('generations'), ylabel('\zeta_0 (rad.)')
% Trip Time
psp4 = subplot(4,1,4);
plot(A(1,:), A(20,:))
title('Trip time')
xlabel('generations'), ylabel('t (s)')

linkaxes([psp1,psp2,psp3,psp4],'x')
linkaxes([sp1,sp2,sp3,sp4,sp5,sp6,psp1,psp2,psp3,psp4],'x')
suptitle('Initial Parameters')

% Final differences, trip time, and annealing graphs
figure(4)
% Final position difference over generations
fsp1 = subplot(2,2,1);
plot(A(1,:), km*A(2,:))
title('Final position difference')
xlabel('generations'), ylabel('distance (km)')
% Final speed difference over generations
fsp2 = subplot(2,2,2);
plot(A(1,:), km*A(3,:))
title('Final speed difference')
xlabel('generations'), ylabel('speed (km/s)')
% Trip time over generations
fsp3 = subplot(2,2,3);
plot(A(1,:), A(20,:))
title('Trip time')
xlabel('generations'), ylabel('t (s)')
% Annealing rate over generations
fsp4 = subplot(2,2,4);
plot(A(1,:), A(19,:))
title('Annealing rate')
xlabel('generations'), ylabel('rate')

linkaxes([sp1,sp2,sp3,sp4,sp5,sp6,psp1,psp2,psp3,psp4,fsp1,fsp2,fsp3,fsp4],'x')


fclose(fin);