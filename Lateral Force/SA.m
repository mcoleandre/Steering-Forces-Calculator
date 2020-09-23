function [ SA_o, SA_i ] = SA(t_f,R_1,l,w_o,R_w,a)
    %Calculating the Steer Angles of the inner and outer wheels
    delta_o = l/(R_1+(t_f/2));
    delta_i = l/(R_1-(t_f/2));
    
    %Calculating the Yaw Rate of the vehicle
    r = (R_w*w_o)/(R_1+(t_f/2));
    
    %Longitudinal Velocity
    v_x = w_o/R_w;

    %Lateral Velocity

    
    %Calculating the Slip Angles of the inner and outer wheels
    SA_o = -delta_o+atan((v_x+r*a)/(v_x-0.5*r*t_f));
    SA_i = -delta_i+atan((v_x+r*a)/(v_x+0.5*r*t_f));
end