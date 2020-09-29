function M_z0 = Self_aligning_moment(gamma, F_z0, V_cx, V_cy)
    
    



    

    D_t0 = F_z * (R_0 / (F_z0 * lambda_Fz0)) * (q_Dz1 + q_Dz2 * d_fz) * lambda_t * sign(V_cx);

    D_t = D_t0 * (1 + q_Dz3 * sin(gamma) + q_Dz4 * (sin(gamma))^2);

    V_c = sqrt((V_cx)^2 + (V_cy)^2);

    alpha_star = (- V_cy) / (abs(V_cx + Epsilon));
    
    d_fz = (F_z - F_z0 * lambda_Fz0) / (F_z0 * lambda_Fz0);
    
    S_Ht = q_Hz1 + (q_Hz2 * d_fz) + (q_Hz3 + q_Hz4 * d_fz)*sin(gamma);
    
    alpha_t = alpha_star + S_Ht;
    
    t_0 = D_t * cos(C_t * atan(B_t * alpha_t - E_t * atan(B_t * alpha_t))) * (V_cx / (V_c + Epsilon_V)); % Pneumatic trail
    
    M_z0_prime = - t_0 * F_y0;
    
    M_z0 = M_z0_prime + M_zr0;
end





