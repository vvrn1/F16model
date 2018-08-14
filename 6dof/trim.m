function [State, Action] = trim(X,Y,velocity,altitude,psi)
close all;
tic
global  State0 dt


%% Trim aircraft to desired altitude and velocity
%%
% altitude = input('Enter the altitude for the simulation (m)  :  ');
% velocity = input('Enter the velocity for the simulation (m/s):  ');
% simtime = input('Enter length of simulation time (s):  ');
% altitude = 6000;
% velocity = 160;
%% Initial Conditions for trim routine.
%% The following values seem to trim to most
%% flight condition.  If the F16 does not trim
%% Change these values.
beta = 0;         % -
elevator = 0*pi/180;       % elevator, rad
alpha = 5*pi/180;           % AOA, rad
rudder = 0;         % rudder angle, rad
aileron = 0;         % aileron, rad
dth = 0.24;
power0 = tgear(dth);
% psi = pi;
phi = 0;
theta = alpha;
X0  = X ;
Y0 = Y ;
Z0 = -altitude;
p_body0 = 0 ;
q_body0 = 0;
r_body0 = 0;
Vt0 = velocity;
% q0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
% q1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);
% q2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
% q3 = -cos(psi/2)*sin(theta/2)*sin(phi/2) + sin(psi/2)*cos(theta/2)*cos(phi/2);

% Initial Guess for free parameters
State0 = [X0, Y0, Z0, Vt0, alpha, beta, phi, theta, psi, p_body0 ,q_body0 ,r_body0, power0];
UX0 = [alpha,beta,aileron,elevator, rudder,dth];
dt =0.1;

% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+04,'MaxIter',1e+04);

% iter = 1;
% cost =1;
% while (iter == 1)
% NextState = F16model(State,Action,dt);
% State = NextState;
% cost = norm(NextState - State);
   
[UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('F16modeltrim',UX0,OPTIONS);
cost = F16modeltrim(UX)
State = [X0, Y0, Z0, Vt0, UX(1), UX(2), cos(UX(1)/2), 0, sin(UX(1)/2), 0, p_body0 ,q_body0 ,r_body0, power0]
Action = [UX(3),UX(4),UX(5),UX(6)]
%     [cost, Xdot, xu, uu] = trim_fun(UX);
    
%     disp('Trim Values and Cost:');
%     disp(['cost   = ' num2str(cost)])
%     disp(['dth    = ' num2str(uu(1)) ' -'])  
%     disp(['elev   = ' num2str(uu(2)*180/pi) ' deg'])
%     disp(['ail    = ' num2str(uu(3)*180/pi) ' deg'])
%     disp(['rud    = ' num2str(uu(4)*180/pi) ' deg'])
%     disp(['alpha  = ' num2str(xu(3)*180/pi) ' deg'])
%     disp(['dLEF   = ' num2str(uu(5)*180/pi) ' deg'])
%     disp(['Vel.   = ' num2str(xu(1)) ' m/s']) 
%     disp(['pow    = ' num2str(xu(14)) ' %']) 
%         iter = 0;
%     feval('F16_trim', [], [], [], 'term');
%     UX0 = UX;
% end

% alpha_trim = xu(3);
% V_trim = xu(1);
% h_trim = altitude;
% pow_trim = xu(14);
% de_trim = uu(2);
% del_trim = uu(2)*180/pi;
% der_trim = uu(2)*180/pi;
% da_trim = uu(3)*180/pi;
% dr_trim = uu(4)*180/pi;
% dlef_trim = uu(5);
% init_x = [xu(1:14);uu(7);uu(8)];
% %init_x = [V_trim 0 alpha_trim, 0 alpha_trim 0, 0 0 0, 0 0 -h_trim, pow_trim];
% pitch_trim = alpha_trim*180/pi;
% yaw_trim =0;
% roll_trim =0;

disp('------------------------------------------');
disp('Ready');
toc