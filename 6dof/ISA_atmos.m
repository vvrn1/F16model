function [mach,qbar,  rho,grav] = ISA_atmos(alt,vt)  
rho0 = 1.225;
Re = 6371000;
R = 287.05;
T0 = 288.15;
g0 = 9.80665;
gamma = 1.4;
if alt >= 11000.0
    temp = 216.65;
elseif alt < 11000.0
    temp = T0 - 0.0065 * alt;
end
rho = rho0 * exp((-g0 / (R * temp)) * alt);
mach = vt / sqrt(gamma * R * temp);
qbar = 0.5 * rho * vt * vt;
grav = g0*(Re*Re/((Re+alt)*(Re+alt)));
