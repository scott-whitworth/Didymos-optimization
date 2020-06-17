function [] = individualProgress(n)
% individualProgress: plots generational data from a binary file
    % Input:
        % n: the number of "individualProgress" binary file to be read in 
    % Output:
        % N/A: see figures

fin = fopen(join(['individualDifference',num2str(n),'.bin']));
A = fread(fin, [Inf, 20], 'double');

au = 1.49587870691e11;
km = au/1e3;

% Position graphs
figure(1)
% R initial over generations
subplot(2,3,1)
plot(A(:,1), km*A(:,10))
title('Radius')
xlabel('generations'), ylabel('r_i (km)')
% Theta initial over generations
subplot(2,3,2)
plot(A(:,1), A(:,11))
title('Direction')
xlabel('generations'), ylabel('\theta_i (rad.)')
% Z initial over generations
subplot(2,3,3)
plot(A(:,1), km*A(:,12))
xlabel('generations'), ylabel('z_i (km)')
title('Elevation')
% R final over generations
subplot(2,3,4)
plot(A(:,1), km*A(:,4))
xlabel('generations'), ylabel('r_f (km)')
% Theta final over generations
subplot(2,3,5)
plot(A(:,1), A(:,5))
xlabel('generations'), ylabel('\theta_f (rad.)')
% Z final over generations
subplot(2,3,6)
plot(A(:,1), km*A(:,6))
xlabel('generations'), ylabel('z_f (km)')

% Velocity graphs
figure(2)
% VR initial over generations
subplot(2,3,1)
plot(A(:,1), km*A(:,13))
title('Radial velocity')
xlabel('generations'), ylabel('v_r_i (km/s)')
% VTh initial over generations
subplot(2,3,2)
plot(A(:,1), km*A(:,14))
title('Angular velocity')
xlabel('generations'), ylabel('rv_\theta_i (km/s)')
% VZ initial over generations
subplot(2,3,3)
plot(A(:,1), km*A(:,15))
title('Upward velocity')
xlabel('generations'), ylabel('v_z_i (km/s)')
% VR final over generations
subplot(2,3,4)
plot(A(:,1), km*A(:,7))
xlabel('generations'), ylabel('v_r_f (km/s)')
% VTh final over generations
subplot(2,3,5)
plot(A(:,1), km*A(:,8))
xlabel('generations'), ylabel('v_\theta_f (km/s)')
% VZ final over generations
subplot(2,3,6)
plot(A(:,1), km*A(:,9))
xlabel('generations'), ylabel('v_z_f (km/s)')

% Initial parameter graphs
figure(3)
% Alpha over generations
subplot(3,1,1)
plot(A(:,1), A(:,16))
title('Initial position angle')
xlabel('generations'), ylabel('\alpha_0 (rad.)')
% Beta over generations
subplot(3,1,2)
plot(A(:,1), A(:,17))
title('Initial in-plane velocity angle')
xlabel('generations'), ylabel('\beta_0 (rad.)')
% Zeta over generations
subplot(3,1,3)
plot(A(:,1), A(:,18))
title('Initial out-of-plane velocity angle')
xlabel('generations'), ylabel('\zeta_0 (rad.)')

% Final differences, trip time, and annealing graphs
figure(4)
% Final position difference over generations
subplot(2,2,1)
plot(A(:,1), km*A(:,2))
title('Final position difference')
xlabel('generations'), ylabel('distance (km)')
% Final speed difference over generations
subplot(2,2,2)
plot(A(:,1), km*A(:,3))
title('Final speed difference')
xlabel('generations'), ylabel('speed (km/s)')
% Trip time over generations
subplot(2,2,3)
plot(A(:,1), A(:,20))
title('Trip time')
xlabel('generations'), ylabel('t (s)')
% Annealing rate over generations
subplot(2,2,4)
plot(A(:,1), A(:,19))
title('Annealing rate')
xlabel('generations'), ylabel('rate')

fclose(fin);
