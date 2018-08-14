function [cost] = F16modeltrim(UX)
%%
global State0 dt
global f16engine
F16parameter;
F16aerodata;
F16enginedata;


%%
%%
X_earth = State0(1);
Y_earth = State0(2);
Z_earth = State0(3);
Vt = State0(4);
phi = State0(7);
theta = UX(1);
psi = State0(9);
alpha = UX(1);
beta = UX(2);
% q0= State0(7);
% q1 = State0(8);
% q2 = State0(9);
% q3 = State0(10);
q0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
q1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);
q2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
q3 = -cos(psi/2)*sin(theta/2)*sin(phi/2) + sin(psi/2)*cos(theta/2)*cos(phi/2);
p_body = State0(10);
q_body = State0(11);
r_body = State0(12);
% da = State(13);
% de = State(13);
% dr = State(13);
% dp = State(13);
% power = State0(14);
da = UX(3);
de = UX(4);
dr = UX(5);
dp = UX(6);
power = tgear(dp);
%%
if beta > 30*pi/180
    beta = 30*pi/180;
elseif beta < -30*pi/180
    beta = -30*pi/180;
end;

% elevator limits
if de > 25*pi/180
    de = 25*pi/180;
elseif de < -25*pi/180
    de = -25*pi/180;
end;

% angle of attack limits
  if alpha > 90*pi/180
    alpha = 90*pi/180;
  elseif alpha < -20*pi/180
    alpha = -20*pi/180;
  end


%  Aileron limits
if da > 21.5*pi/180
    da = 21.5*pi/180;
elseif da < -21.5*pi/180
    da = -21.5*pi/180;
end;

% Rudder limits
if dr > 30*pi/180
    dr = 30*pi/180;
elseif dr< -30*pi/180
    dr = -30*pi/180;
end;

   

% dth limits
if dp > 1
    dp = 1;
elseif dp < 0
    dp = 0;
end;

% Verify that the calculated leading edge flap
% have not been violated.

%%
[mach,qbar, rho, g] = ISA_atmos(-Z_earth,Vt) ;
%%
pdot = power_dot(dp,power);
thrust = engine_power(mach,-Z_earth,power);
%%  Body velocity components 
u_body = Vt * cos(alpha) * cos(beta);
v_body = Vt * sin(beta);
w_body = Vt * sin(alpha) * cos(beta);

%% Transformation rad to deg for the lookup tables 
ALPHA = alpha * rtd;
BETA = beta * rtd;
DE = de * rtd;
DA = da * rtd;
DR = dr * rtd;
DLEF = 1.38 *ALPHA -9.05 * qbar /(101325*rho/1.225) +1.45;
% DLEF = dlef * rtd;
if (DLEF > 25*pi/180)
    DLEF = 25*pi/180;
elseif (DLEF < 0)
    DLEF = 0;
end;
%% LEF tables are only valid up until alpha = 45 degrees
if (ALPHA > 45)
    alpha_lef = 45;
else
    alpha_lef = ALPHA;
end
%%  Limit of elevator deflection is 25 degrees
if (DE > 25)
    DE = 25;
elseif (DE < -25)
    DE = -25;
end
    
%% Normalizing the control deflections 
da_norm = DA/21.5;
dr_norm = DR/30.0;
dlef_norm = (1 - DLEF/25.0);
%% hifi lookup tables 
% CX
CX0 = interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.CX,BETA,ALPHA,DE);
delta_CX_lef = interp2(f16data.beta,f16data.alpha2,f16data.CX_lef,BETA,alpha_lef)-interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.CX,BETA,ALPHA,0);
CXq = interp1(f16data.alpha1,f16data.CXq,ALPHA);
dCXq_lef = interp1(f16data.alpha2,f16data.dCXq_lef,alpha_lef);
% CZ
CZ0 = interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.CZ,BETA,ALPHA,DE);

delta_CZ_lef = interp2(f16data.beta,f16data.alpha2,f16data.CZ_lef,BETA,alpha_lef) - interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.CZ,BETA,ALPHA,0);

CZq = interp1(f16data.alpha1,f16data.CZq,ALPHA);

dCZq_lef = interp1(f16data.alpha2,f16data.dCZq_lef,alpha_lef);
% Cm
Cm0 = interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.Cm,BETA,ALPHA,DE);

delta_Cm_lef = interp2(f16data.beta,f16data.alpha2,f16data.Cm_lef,BETA,alpha_lef)- ...
    interp3(f16data.beta,f16data.alpha1,f16data.de1,f16data.Cm,BETA,ALPHA,0);

Cmq = interp1(f16data.alpha1,f16data.Cmq,ALPHA);

dCmq_lef = interp1(f16data.alpha2,f16data.dCmq_lef,alpha_lef);

dCm = interp1(f16data.alpha1,f16data.dCm,ALPHA);

dCm_ds = interp2(f16data.de3,f16data.alpha1,f16data.dCm_ds,DE,ALPHA);
% CY
CY0 = interp2(f16data.beta,f16data.alpha1,f16data.CY,BETA,ALPHA);

delta_CY_lef= interp2(f16data.beta,f16data.alpha2,f16data.CY_lef,BETA,alpha_lef)-CY0;

delta_CY_da20 = interp2(f16data.beta,f16data.alpha1,f16data.CY_da20,BETA,ALPHA)-CY0;

delta_CY_da20lef = interp2(f16data.beta,f16data.alpha2,f16data.CY_da20lef,BETA,alpha_lef)- ...
    interp2(f16data.beta,f16data.alpha2,f16data.CY_lef,BETA,alpha_lef) - ...
   delta_CY_da20 ;

delta_CY_dr30= interp2(f16data.beta,f16data.alpha1,f16data.CY_dr30,BETA,ALPHA) - CY0;

CYr = interp1(f16data.alpha1,f16data.CYr,ALPHA);

dCYr_lef = interp1(f16data.alpha2,f16data.dCYr_lef,alpha_lef);

CYp = interp1(f16data.alpha1,f16data.CYp,ALPHA);

dCYp_lef = interp1(f16data.alpha2,f16data.dCYp_lef,alpha_lef);

% Cn
Cn0 = interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cn,BETA,ALPHA,DE);

delta_Cn_lef = interp2(f16data.beta,f16data.alpha2,f16data.Cn_lef,BETA,alpha_lef) - ...
   interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cn,BETA,ALPHA,0) ;

delta_Cn_da20 = interp2(f16data.beta,f16data.alpha1,f16data.Cn_da20,BETA,ALPHA) - ...
    interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cn,BETA,ALPHA,0);

delta_Cn_da20lef = interp2(f16data.beta,f16data.alpha2,f16data.Cn_da20lef,BETA,alpha_lef) - ...
     interp2(f16data.beta,f16data.alpha2,f16data.Cn_lef,BETA,alpha_lef) - ...
    delta_Cn_da20  ;

delta_Cn_dr30 = interp2(f16data.beta,f16data.alpha1,f16data.Cn_dr30,BETA,ALPHA) - ...
    interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cn,BETA,ALPHA,0) ;

Cnr = interp1(f16data.alpha1,f16data.Cnr,ALPHA);

dCnbeta = interp1(f16data.alpha1,f16data.dCnbeta,ALPHA);

dCnr_lef = interp1(f16data.alpha2,f16data.dCnr_lef,alpha_lef);

Cnp = interp1(f16data.alpha1,f16data.Cnp,ALPHA);

dCnp_lef = interp1(f16data.alpha2,f16data.dCnp_lef,alpha_lef);
% Cl
Cl0 = interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cl,BETA,ALPHA,DE);

delta_Cl_lef = interp2(f16data.beta,f16data.alpha2,f16data.Cl_lef,BETA,alpha_lef) - ...
    interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cl,BETA,ALPHA,0);

delta_Cl_da20 = interp2(f16data.beta,f16data.alpha1,f16data.Cl_da20,BETA,ALPHA) - ...
    interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cl,BETA,ALPHA,0);

delta_Cl_da20lef = interp2(f16data.beta,f16data.alpha2,f16data.Cl_da20lef,BETA,alpha_lef)- ...
   interp2(f16data.beta,f16data.alpha2,f16data.Cl_lef,BETA,alpha_lef) - ...
   delta_Cl_da20 ;

delta_Cl_dr30 = interp2(f16data.beta,f16data.alpha1,f16data.Cl_dr30,BETA,ALPHA) - ...
    interp3(f16data.beta,f16data.alpha1,f16data.de2,f16data.Cl,BETA,ALPHA,0);

Clr = interp1(f16data.alpha1,f16data.Clr,ALPHA);

dClbeta = interp1(f16data.alpha1,f16data.dClbeta,ALPHA);

dClr_lef = interp1(f16data.alpha2,f16data.dClr_lef,alpha_lef);

Clp = interp1(f16data.alpha1,f16data.Clp,ALPHA);

dClp_lef = interp1(f16data.alpha2,f16data.dClp_lef,alpha_lef);
%%   Total force coefficients 
%  Cx_tot 
CX_tot = CX0 + delta_CX_lef  * dlef_norm + ...
    (F16data.cref/(2*Vt))*(CXq + dCXq_lef * dlef_norm) * q_body;
% Cy_tot
CY_tot = CY0 + delta_CY_lef * dlef_norm+ ...
    (delta_CY_da20 + delta_CY_da20lef  * dlef_norm) * da_norm+ ...
    delta_CY_dr30 * dr_norm+ ...
    (F16data.bref / (2*Vt))*(CYr + dCYr_lef * dlef_norm) * r_body+ ...
    (F16data.bref/(2*Vt))*(CYp + dCYp_lef * dlef_norm) * p_body;
%     Cz_tot 
CZ_tot = CZ0 + delta_CZ_lef * dlef_norm + ...
      (F16data.cref/(2*Vt))*(CZq + dCZq_lef * dlef_norm) * q_body;
 %  Cl_tot 
 Cl_tot = Cl0 + delta_Cl_lef * dlef_norm ... 
     + (delta_Cl_da20 + delta_Cl_da20lef * dlef_norm) * da_norm ...
     + delta_Cl_dr30 * dr_norm ...
     + (F16data.bref / ((2*Vt))*(Clr + dClr_lef * dlef_norm)) * r_body ...
     + ((F16data.bref / (2*Vt)) * (Clp + dClp_lef * dlef_norm)) * p_body ...
     + dClbeta* BETA;
  %     Cm_tot
  Cm_tot = Cm0 * 1.0 + CZ_tot * (F16data.xcgr - F16data.xcg) + delta_Cm_lef * dlef_norm ...
     + (F16data.cref / (2*Vt))*(Cmq + dCmq_lef  * dlef_norm) * q_body ...
     + dCm + dCm_ds ;
 % Cn_tot
     Cn_tot = Cn0 + delta_Cn_lef * dlef_norm ...
     - CY_tot * (F16data.xcgr - F16data.xcg)*(F16data.cref/F16data.bref) ...
     + (delta_Cn_da20 + delta_Cn_da20lef * dlef_norm) * da_norm ...
     + ((F16data.bref / (2*Vt)) * (Cnr + dCnr_lef * dlef_norm))* r_body ...
     + ((F16data.bref / (2*Vt)) * (Cnp + dCnp_lef * dlef_norm)) * p_body ...
     + delta_Cn_dr30  * dr_norm + dCnbeta * beta * rtd;
 %% 	 Total forces 
 Xbar = qbar * F16data.Sref * CX_tot;
 Ybar = qbar * F16data.Sref * CY_tot;
 Zbar = qbar * F16data.Sref * CZ_tot;
 %%    Total moments 
 Lbar = Cl_tot * qbar * F16data.Sref * F16data.bref;
 Mbar = Cm_tot * qbar * F16data.Sref * F16data.cref;
 Nbar = Cn_tot * qbar * F16data.Sref * F16data.bref;
%%   Derivatives 
% q0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
% q1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);
% q2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
% q3 = -cos(psi/2)*sin(theta/2)*sin(phi/2) + sin(psi/2)*cos(theta/2)*cos(phi/2);

u_body_dot = r_body * v_body - q_body * w_body ...
    + (Xbar + thrust) / F16data.mass + 2*(q1*q3 - q0*q2)*g;
v_body_dot = p_body * w_body - r_body * u_body ...
    + Ybar / F16data.mass + 2*(q2*q3 + q0*q1)*g;
w_body_dot = q_body * u_body - p_body * v_body ...
    + Zbar / F16data.mass + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*g;

Vt_dot      = (u_body * u_body_dot + v_body * v_body_dot ...
    + w_body * w_body_dot) / Vt;
beta_dot    = (v_body_dot * Vt - v_body * Vt_dot) ...
    / (Vt * Vt * cos(beta));
alpha_dot   = (u_body * w_body_dot - w_body * u_body_dot) ...
    / (u_body * u_body + w_body * w_body);

q0_dot      = 0.5 * (-p_body * q1 - q_body * q2 - r_body * q3);
q1_dot      = 0.5 * ( p_body * q0 + r_body * q2 - q_body * q3);
q2_dot      = 0.5 * ( q_body * q0 - r_body * q1 + p_body * q3);
q3_dot      = 0.5 * ( r_body * q0 + q_body * q1 - p_body * q2);

%  correction term 
dq = q0 * q0_dot + q1 * q1_dot + q2 * q2_dot + q3 * q3_dot;

q0_dot = q0_dot - dq * q0;
q1_dot = q1_dot - dq * q1;
q2_dot = q2_dot - dq * q2;
q3_dot = q3_dot - dq * q3;

p_body_dot  = (F16data.C1 * r_body + F16data.C2 * p_body) * q_body ...
    +F16data.C3 * Lbar +F16data.C4 * (Nbar + q_body * F16data.heng);
q_body_dot  = F16data.C5 * p_body * r_body ...
    - F16data.C6 * (p_body * p_body - r_body * r_body) ...
    + F16data.C7 * (Mbar - F16data.heng * r_body);
r_body_dot  = (F16data.C8 * p_body - F16data.C2 * r_body) * q_body ...
    + F16data.C4 * Lbar + F16data.C9 * (Nbar + q_body * F16data.heng);

x_earth_dot = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * u_body ...
    + 2*(q1*q2 - q0*q3) * v_body ...
    + 2*(q1*q3 + q0*q2) * w_body;
y_earth_dot = 2*(q1*q2 + q0*q3) * u_body ...
    + (q0*q0 - q1*q1 + q2*q2 - q3*q3) * v_body ...
    + 2*(q2*q3 - q0*q1) * w_body;
z_earth_dot = 2*(q1*q3 - q0*q2) * u_body ...
    + 2*(q2*q3 + q0*q1) * v_body ...
    + (q0*q0 - q1*q1 - q2*q2 + q3*q3) * w_body;
%% intergrate
X_earth = X_earth +x_earth_dot *dt;
Y_earth = Y_earth +y_earth_dot *dt;
Z_earth = Z_earth +z_earth_dot *dt;
Vt = Vt + Vt_dot *dt;
alpha = alpha + alpha_dot*dt;
beta = beta+beta_dot*dt;
q0 = q0 +q0_dot*dt;
q1 = q1 +q1_dot*dt;
q2 = q2 +q2_dot*dt;
q3 = q3 +q3_dot*dt;
p_body = p_body+p_body_dot*dt;
q_body = q_body+q_body_dot*dt;
r_body = r_body+r_body_dot*dt;
power  = power +pdot*dt;
%% normal accelerations 
% ny = Ybar/mass/g;
% nz = -Zbar/mass/g;
%%
NextState(1)= X_earth ;
NextState(2) = Y_earth   ;
NextState(3) = Z_earth;
NextState(4) = Vt;
NextState(5) = alpha;
NextState(6) = beta;
NextState(7) = q0;
NextState(8) = q1;
NextState(9) = q2;
NextState(10) = q3;
NextState(11) = p_body;
NextState(12) = q_body;
NextState(13) = r_body;
% da = State(13);
% de = State(13);
% dr = State(13);
% dp = State(13);
NextState(14) = power;
%%
weight = [  2            ...%Vt_dot
            10           ...%beta_dot
            100           ...%alpha_dot
            10           ...%q0_dot
            10           ...%q1_dot
            10           ...%q2_dot
            10           ...%q3_dot
            10           ...%p_dot
            10           ...%q_dot
            10           ...%r_dot
            0            ...%x_dot
            0            ...%y_dot
            5            ...%z_dot
            50           ...
            ];            %pow_dot
Xdot = [Vt_dot; beta_dot; alpha_dot ; q0_dot;q1_dot;q2_dot;q3_dot;p_body_dot;q_body_dot;r_body_dot;x_earth_dot;y_earth_dot;z_earth_dot;pdot];
cost = weight*(Xdot.*Xdot);
end