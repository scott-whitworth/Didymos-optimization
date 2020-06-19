function [] = thrustProgress(str)
% thrustProgress: plots generational thruster data from a binary file
    % Input:
        % str: the specific generation to be plotted ('Best' or 'Worst')
    % Output:
        % N/A: see figures

fin = fopen(join([str,'ThrustGens.bin']));
A = fread(fin, [16 Inf], 'double');

% Gamma
figure(1)
g0 = subplot(2,4,1);
plot(A(1,:),A(2,:))
xlabel('generations'), ylabel('a_0')
g1 = subplot(2,4,2);
plot(A(1,:),A(3,:))
xlabel('generations'), ylabel('a_1')
g2 = subplot(2,4,3);
plot(A(1,:),A(4,:))
xlabel('generations'), ylabel('b_1')
g3 = subplot(2,4,4);
plot(A(1,:),A(5,:))
xlabel('generations'), ylabel('a_2')
g4 = subplot(2,4,5);
plot(A(1,:),A(6,:))
xlabel('generations'), ylabel('b_2')
g5 = subplot(2,4,6);
plot(A(1,:),A(7,:))
xlabel('generations'), ylabel('a_3')
g6 = subplot(2,4,7);
plot(A(1,:),A(8,:))
xlabel('generations'), ylabel('b_3')
linkaxes([g0 g1 g2 g3 g4 g5 g6], 'x')
suptitle('\gamma coefficients')

% Tau
figure(2)
t0 = subplot(2,2,1);
plot(A(1,:),A(9,:))
xlabel('generations'), ylabel('a_0')
t1 = subplot(2,2,2);
plot(A(1,:),A(10,:))
xlabel('generations'), ylabel('a_1')
t2 = subplot(2,2,3);
plot(A(1,:),A(11,:))
xlabel('generations'), ylabel('b_1')
linkaxes([t0 t1 t2], 'x')
suptitle('\tau coefficients')

% Coast
figure(3)
title('Coast over generations')
c0 = subplot(2,3,1);
plot(A(1,:),A(12,:))
xlabel('generations'), ylabel('a_0')
c1 = subplot(2,3,2);
plot(A(1,:),A(13,:))
xlabel('generations'), ylabel('a_1')
c2 = subplot(2,3,3);
plot(A(1,:),A(14,:))
xlabel('generations'), ylabel('b_1')
c3 = subplot(2,3,4);
plot(A(1,:),A(15,:))
xlabel('generations'), ylabel('a_2')
c4 = subplot(2,3,5);
plot(A(1,:),A(16,:))
xlabel('generations'), ylabel('b_2')
linkaxes([c0 c1 c2 c3 c4], 'x')
suptitle('coast coefficients')

fclose(fin);