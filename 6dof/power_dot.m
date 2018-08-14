function pdot = power_dot(dp,Pa)
if ( dp <= 0.77 )
    Pc = 64.94 * dp;
else
    Pc = 217.38 * dp - 117.38;
end
if (Pc-Pa)<=25
    teng_star = 1.0;
elseif (Pc-Pa)>=50
    teng_star = 0.1;
elseif (Pc-Pa)>25 && (Pc-Pa)<50.0
    teng_star = 1.9 - 0.036*(Pc-Pa);
end
if Pc >=50.0
    if Pa>=50.0
        Pc = Pc;
        teng = 5.0;
    elseif Pa<50
        Pc = 60;
        teng = teng_star;
    end
elseif Pc<50.0
    if Pa>=50.0
        Pc = 40;
        teng = 5.0;
    elseif Pa<50
        Pc = Pc;
        teng = teng_star;
    end
end
pdot = teng * (Pc-Pa);

    
