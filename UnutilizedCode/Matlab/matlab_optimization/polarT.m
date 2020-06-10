function [R,V] = polarT(x,y,z,vx,vy,vz)

[theta, rho, z] = cart2pol(x,y,z);
R=[rho,theta,z];
vv=[vx;vy;vz];
% Transformation to cylindrical for asteroid
vT=[cos(theta)  sin(theta) 0
         -sin(theta) cos(theta) 0
         0           0          1];
% Velocities in cylindrical coordinates for earth
V=vT*vv;
V=transpose(V);
end

