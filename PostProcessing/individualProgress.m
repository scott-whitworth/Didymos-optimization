function [] = individualProgress(filename, n)
% individualProgress: plots generational data from a binary file
    % Input:
        % filename: the name of the data file (ending in .bin)
        % n: the number of rows to read and plot
    % Output:
        % N/A: see figures

fin = fopen(filename);
[gen, dp, dv, rf, thf, zf, vrf, vthf, vzf, ri, thi, zi, vri, vthi, vzi, a, b, z, an, t] = fread(fin, [20, n], 'double');

% Position graphs
figure(1)
% posDiff over generations
subplot(2,2,1)
plot(gen, dp)
title('Final position difference')
xlabel('generations'), ylabel('distance (a.u.)')
% R (initial and final) over generations
subplot(2,2,2)
plot(gen, ri, gen, rf)
legend('initial', 'final')
title('Radius')
xlabel('generations'), ylabel('r (a.u.)')
% Theta (initial and final) over generations
subplot(2,2,3)
plot(gen, thi, gen, thf)
title('Direction')
xlabel('generations'), ylabel('\theta (rad.)')
% Z (initial and final) over generations
subplot(2,2,4)
plot(gen, zi, gen, zf)
title('Elevation')
xlabel('generations'), ylabel('z (a.u.)')

% Velocity graphs
figure(2)
% velDiff over generations
subplot(2,2,1)
plot(gen, dv)
title('Final speed difference')
xlabel('generations'), ylabel('speed (a.u./s)')
% VR (initial and final) over generations
subplot(2,2,2)
plot(gen, vri, gen, vrf)
legend('initial', 'final')
title('Radial velocity')
xlabel('generations'), ylabel('v_r (a.u./s)')
% VTh (initial and final) over generations
subplot(2,2,3)
plot(gen, vthi, gen, vthf)
title('Angular velocity')
xlabel('generations'), ylabel('rv_\theta (a.u./s)')
% VZ (initial and final) over generations
subplot(2,2,4)
plot(gen, vzi, gen, vzf)
title('Upward velocity')
xlabel('generations'), ylabel('v_z (a.u./s)')

% Initial parameter graphs
figure(3)
% Alpha (initial) over generations
subplot(3,1,1)
plot(gen, a)
title('Initial alpha')
xlabel('generations'), ylabel('\alpha_0 (rad.)')
% Beta (initial) over generations
subplot(3,1,2)
plot(gen, b)
title('Initial beta')
xlabel('generations'), ylabel('\beta_0 (rad.)')
% Zeta (initial) over generations
subplot(3,1,3)
plot(gen, z)
title('Initial zeta')
xlabel('generations'), ylabel('\zeta_0 (rad.)')

% Trip time and annealing graphs
figure(4)
% Trip time over generations
subplot(2,1,1)
plot(gen, t)
title('Trip time')
xlabel('generations'), ylabel('t (s)')
% Annealing rate over generations
subplot(2,1,2)
plot(gen, an)
title('Annealing rate')
xlabel('generations'), ylabel('rate')

fclose(fin);