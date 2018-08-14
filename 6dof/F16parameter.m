F16data.mass = 9295.44 ;
F16data.Ixx= 12874.8;
F16data.Iyy   =      75673.6;
F16data.Izz      =   85552.1;
F16data.Ixz        = 1331.4;
F16data.Sref    =    27.87;
F16data.bref     =   9.144;
F16data.cref     =   3.45;
F16data.xcg    =     0.3;
F16data. xcgr     =   0.35;
F16data. heng     =   216.9 ;


rtd      =     57.29577951;
dtr       =    0.017453293;
Pi			=  3.141592654;
F16data.Gamma = F16data.Ixx * F16data.Izz - (F16data.Ixz  * F16data.Ixz);
F16data.C1 = ((F16data.Iyy - F16data.Izz)  * F16data.Izz  - (F16data.Ixz * F16data.Ixz))/ F16data.Gamma;
F16data.C2 = ((F16data.Ixx - F16data.Iyy + F16data.Izz ) * F16data.Ixz ) / F16data.Gamma;
F16data.C3 = F16data.Izz / F16data.Gamma;
F16data.C4 = F16data.Ixz / F16data.Gamma;
F16data.C5 = (F16data.Izz - F16data.Ixx) / F16data.Iyy;
F16data.C6 = F16data.Ixz / F16data.Iyy ;
F16data.C7 = 1 / F16data.Iyy;
F16data.C8 = (F16data.Ixx * (F16data.Ixx - F16data.Iyy ) + F16data.Ixz * F16data.Ixz) /F16data. Gamma;
F16data.C9 = F16data.Ixx / F16data.Gamma;