function [ Dm_fr, Dm_fl, Dm_rr, Dm_rl, DF_zfr, DF_zfl, DF_zrr, DF_zrl ] = DM(m,WT_y,p,g)
    %Static Mass
    sm_fr = (m/2)*(p/100);
    sm_fl = (m/2)*(p/100);
    sm_rr = (m/2)*(1-p/100);
    sm_rl = (m/2)*(1-p/100);

    %Dynamic Mass
    Dm_fr = sm_fr + WT_y*(p/100) - (WT_xi)/2;
    Dm_fl = sm_fl - WT_y*(p/100) - (WT_xi)/2;
    Dm_rr = sm_rr + WT_y*(p/100) + (WT_xi)/2;
    Dm_rl = sm_rl - WT_y*(p/100) + (WT_xi)/2;
    
    %Dynamic Normal Loads
    DF_zfr = Dm_fr*g;
    DF_zfl = Dm_fl*g;
    DF_zrr = Dm_rr*g;
    DF_zrl = Dm_rl*g;
end