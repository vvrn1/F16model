function thrust = engine_power(Ma,alt,Pa)
global f16engine
alt = alt / 0.3048;
T_Idle = interp2(f16engine.alt,f16engine.Ma,f16engine.Idle,alt,Ma);
T_Mil = interp2(f16engine.alt,f16engine.Ma,f16engine.Mil,alt,Ma);
T_Max = interp2(f16engine.alt,f16engine.Ma,f16engine.Max,alt,Ma);
if Pa<50
    thrust = T_Idle +(T_Mil - T_Idle)*Pa/50;
elseif Pa>=50
    thrust = T_Mil +(T_Max - T_Mil)*(Pa-50)/50;
end
thrust = thrust*4.4482216;
    
    
