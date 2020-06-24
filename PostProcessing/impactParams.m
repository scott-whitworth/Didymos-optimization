% impactParams(): script to calculate position and velocity vector
% components in cylindrical coordinates at impact using Cartesian vector
% components retrieved from https://ssd.jpl.nasa.gov/sbdb.cgi 
clc,clear
day2sec = 3600*24;

% Asteroid position
x_a = 1.022170298157143; y_a = 1.639822634888091e-1;
z_a = -5.541927402432133e-2;
[th_a, r_a] = cart2pol(x_a,y_a);
r_a = vpa(r_a), th_a = vpa(th_a), z_a = vpa(z_a)
% Asteroid velocity
vx_a = -5.273553581680452e-3/day2sec; vy_a = 1.904434089383773e-2/day2sec;
vz_a = 6.305942169924792e-4/day2sec;
vr_a = sqrt((cos(th_a)*vx_a)^2 + (sin(th_a)*vy_a)^2); 
vth_a = sqrt((-sin(th_a)*vx_a)^2 + (cos(th_a)*vy_a)^2); 
vr_a = vpa(vr_a), vth_a = vpa(vth_a), vz_a = vpa(vz_a)

% Earth position
x_e = 9.932443486994552e-1; y_e = 1.277210104878882e-1;
z_e = -1.077701132014302e-5;
[th_e, r_e] = cart2pol(x_e,y_e);
r_e = vpa(r_e), th_e = vpa(th_e), z_e = vpa(z_e)
% Earth velocity
vx_e = -2.480282587602037e-3/day2sec; vy_e = 1.700148842113997e-2/day2sec;
vz_e = -2.573773002330340e-7/day2sec;
vr_e = sqrt((cos(th_e)*vx_e)^2 + (sin(th_e)*vy_e)^2); 
vth_e = sqrt((-sin(th_e)*vx_e)^2 + (cos(th_e)*vy_e)^2);
vr_e = vpa(vr_e), vth_e = vpa(vth_e), vz_e = vpa(vz_e)