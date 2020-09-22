%Right Hand Turn
%Outside Wheel is Left
%Inside Wheel is Right
clc;
clear;
close all;
check = false;

%Inputs
t_f
R_1
l
a
w_o
w_i
R_w
m
p
CG
a_yi
a_xi
g
e
gamma

%Pacejka 2002 Lateral Force Model Parameters
p_Cy1
p_Dy1
p_Dy2
p_Dy3
p_Ey1
p_Ey2
p_Ey3
p_Ey4
p_Ky1
p_Ky2
p_Ky3
p_Hy1
p_Hy2
p_Hy3
p_Vy1
p_Vy2
p_Vy3
p_Vy4

%Pacejka 2002 Self-Aligning Moment Model Parameters
q_By1
q_By2
q_By3
q_By5
q_By6
q_By9
q_By10
q_Cy1
q_Dy1
q_Dy2
q_Dy3
q_Dy4
q_Dy6
q_Dy7
q_Dy8
q_Dy9
q_Ey1
q_Ey2
q_Ey3
q_Ey4
q_Ey5
q_Hy1
q_Hy2
q_Hy3
q_Hy4

%Scaling Factors
lambda_Fz0
lambda_muxy
lambda_muv
lambda_Cy
lambda_Ey
lambda_Kya
lambda_Kygamma
lambda_Mr
lambda_Kzy
lambda_Kz

%Steer Angles
delta_o = l/(R_1+(t_f/2));
delta_i = l/(R_1-(t_f/2));

%Yaw Rate
r = (R_w*w_o)/(R_1+(t_f/2));

%Longitudinal Velocity
v_x = w_o/R_w;

%Lateral Velocity


%Initial Weight Transfer
WT_yi = 0;
WT_xi = 0;
WT_y = WT_yi;

%Static Mass
sm_fr = (m/2)*(p/100);
sm_fl = (m/2)*(p/100);
sm_rr = (m/2)*(1-p/100);
sm_rl = (m/2)*(1-p/100);

%Static Normal Load
sF_zfr = sm_fr*g;
sF_zfl = sm_fl*g;
sF_zrr = sm_rr*g;
sF_zrl = sm_rl*g;

%Dynamic Mass
Dm_fr = sm_fr + WT_yi*(p/100) - (WT_xi)/2;
Dm_fl = sm_fl - WT_yi*(p/100) - (WT_xi)/2;
Dm_rr = sm_rr + WT_yi*(p/100) + (WT_xi)/2;
Dm_rl = sm_rl - WT_yi*(p/100) + (WT_xi)/2;

%Dynamic Normal Loads
DF_zfr = Dm_fr*g;
DF_zfl = Dm_fl*g;
DF_zrr = Dm_rr*g;
DF_zrl = Dm_rl*g;

%General Equations
%Tangent of slip angle (Lateral slip)
alpha = -v_y/abs(v_x+e);

%Digressive Friction Factor
lambda_muV = (A_mu*lambda_muxy)/(1+(A_mu-1)*lambda_muxy);

while check == false
    %Dynamic Mass
    Dm_fr = sm_fr + WT_yi*(p/100) - (WT_xi)/2;
    Dm_fl = sm_fl - WT_yi*(p/100) - (WT_xi)/2;
    Dm_rr = sm_rr + WT_yi*(p/100) + (WT_xi)/2;
    Dm_rl = sm_rl - WT_yi*(p/100) + (WT_xi)/2;
    
    %Dynamic Normal Loads
    DF_zfr = Dm_fr*g;
    DF_zfl = Dm_fl*g;
    DF_zrr = Dm_rr*g;
    DF_zrl = Dm_rl*g;
    
    %Normalized change in vertical load
    df_zfr = (DF_zfr-DF_zfr*lambda_Fz0)/(DF_zfr*lambda_Fz0);
    df_zfl = (DF_zfl-DF_zfl*lambda_Fz0)/(DF_zfl*lambda_Fz0);

    %Inner Wheel
    C_yfr = p_cy1*lambda_Cy;
    D_yfr = mu_yfr*DF_zfr;
    mu_yfr = (p_Dy1 + p_Dy2*df_zfr)*(1-(p_Dy3*(sin(gamma))^2)*lambda_muV);
    K_yalpha0fr = p_Ky1*DF_zfr*lambda_Fz0*sin(2*atan(DF_zfr/(p_Ky2*DF_zfr*lambda_Fz0)))*lambda_Ey;
    K_yalphafr = K_yalpha0fr*(1-p_Ky3*(sin(gamma))^2);
    B_yfr = K_yalphafr/(C_yfr*D_yfr+e);
    S_Hyfr = (p_Hy1 + p_Hy2*df_zfr)*lambda_muxy+p_Hy3*sin(gamma)*lambda_Kygamma;
    alpha_yfr = alpha + S_Hyfr;
    S_Vyfr = DF_zfr*((p_Vy1+p_Vy2*df_zfr)*lambda_muxy+(p_Vy3+p_Vy4*df_zfr)*sin(gamma)*lambda_Kygamma)*lambda_muV;
    E_yfr = (p_Ey1 + p_Ey2*df_zfr)*(1-(p_Ey3 + p_Ey4*sin(gamma))*sign(alpha_yfr))*lambda_Ey;
    K_ygamma0fr = (p_Hy3*K_yalpha0fr+DF_zfr*(p_V3_p_Vy3*df_zfr))*lambda_Kygamma;
    
    F_y0fr = D_yfr*sin(C_yfr*atan(B_yfr*alpha_yfr-E_yfr(B_yfr*alpha_yfr-atan(B_yfr*alpha_yfr))))+S_vyfr;
    
    %Outer Wheel
    C_yfl = p_cy1*lambda_Cy;
    D_yfl = mu_yfr*DF_zfl;
    mu_yfl = (p_Dy1 + p_Dy2*df_zfl)*(1-(p_Dy3*(sin(gamma))^2)*lambda_muV);
    K_yalpha0fl = p_Ky1*DF_zfl*lambda_Fz0*sin(2*atan(DF_zfl/(p_Ky2*DF_zfl*lambda_Fz0)))*lambda_Ey;
    K_yalphafl = K_yalpha0fl*(1-p_Ky3*(sin(gamma))^2);
    B_yfl = K_yalphafl/(C_yfl*D_yfl+e);
    S_Hyfl = (p_Hy1 + p_Hy2*df_zfl)*lambda_muxy+p_Hy3*sin(gamma)*lambda_Kygamma;
    alpha_yfl = alpha + S_Hyfl;
    S_Vyfl = DF_zfl*((p_Vy1+p_Vy2*df_zfl)*lambda_muxy+(p_Vy3+p_Vy4*df_zfl)*sin(gamma)*lambda_Kygamma)*lambda_muV;
    E_yfl = (p_Ey1 + p_Ey2*df_zfl)*(1-(p_Ey3 + p_Ey4*sin(gamma))*sign(alpha_y))*lambda_Ey;
    K_ygamma0fl = (p_Hy3*K_yalpha0fl+DF_zfl*(p_V3_p_Vy3*df_zfl))*lambda_Kygamma;
    
    F_y0fl = D_yfl*sin(C_yfl*atan(B_yfl*alpha_yfl-E_yfl(B_yfl*alpha_yfl-atan(B_yfl*alpha_yfl))))+S_vyfl;
    
    %Lateral Acceleration
    a_yfr = F_y0fr/Dm_fr;
    a_yfl = F_y0fl/Dm_fl;
    a_y = (a_yfr+a_yfl)/2;
    
    %New Weight Transfer
    WT_ynew = (m*CG*a_y)/t_f;
    
    if (WT_y-0.005) < WT_ynew && WT_ynew > (WT_y+0.005)
       check = true;
    end    
end