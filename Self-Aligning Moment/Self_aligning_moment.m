function [M_z0, K_zalpha0, K_zgamma0] = Self_aligning_moment(gamma, F_z0, F_y0, V_cx, V_cy)
    
    d_fz = (F_z - F_z0 * lambda_Fz0) / (F_z0 * lambda_Fz0);
    
    V_c = sqrt((V_cx)^2 + (V_cy)^2);
    
    D_r = F_z * R_0 * ((q_Dz6 + q_Dz7 * d_fz) * lambda_Mr + (q_Dz8 + q_Dz9 * d_fz) * sin(gamma) * lambda_Kzy) * (V_cx / (V_c + Epsilon_V)) * lambda_muy * sign(V_cx);

    C_r = 1;

    B_r = q_Bz10 * B_y *C_y;

    K_yalpha_prime = K_yalpha + Epsilon_K;

    S_Hf = S_Hy + (S_Vy / K_yalpha_prime);
    
    S_Ht = q_Hz1 + (q_Hz2 * d_fz) + (q_Hz3 + q_Hz4 * d_fz)*sin(gamma);
    
    alpha_star = (- V_cy) / (abs(V_cx + Epsilon));

    alpha_r = alpha_star + S_Hf;

    alpha_t = alpha_star + S_Ht;
    
    C_t = q_Cz1;

    B_t = (q_Bz1 + q_Bz2 * d_fz + q_Bz3 * (d_fz)^2) * (1 + q_Bz4 * sin(gamma) + q_Bz5 * abs(sin(gamma))) * lambda_Kyalpha * lambda_muy;

    E_t = (q_Ez1 + q_Ez2 * d_fz + q_Ez3 * (d_fz)^2) * (1 + (q_Ez4 + q_Ez5 * sin(gamma)) * (2 / pi ) * atan(B_t * C_t * alpha_t));

    D_t0 = F_z * (R_0 / (F_z0 * lambda_Fz0)) * (q_Dz1 + q_Dz2 * d_fz) * lambda_t * sign(V_cx);

    D_t = D_t0 * (1 + q_Dz3 * sin(gamma) + q_Dz4 * (sin(gamma))^2);
    
    t_0 = D_t * cos(C_t * atan(B_t * alpha_t - E_t * atan(B_t * alpha_t))) * (V_cx / (V_c + Epsilon_V)); % Pneumatic trail
    
    M_z0_prime = - t_0 * F_y0;
    
    M_zr0 = D_r * cos(C_r * atan(B_r * alpha_r)); % This is the residual self aligning moment

    
    M_z0 = M_z0_prime + M_zr0; % This is the self aligning moment
    K_zalpha0 = D_t0 * K_yalpha0; % This is the slip angle aligning stiffness
    K_zgamma0 = F_z * R_0 * (q_Dz8 + q_Dz9 * d_fz) * lambda_Kzgamma - D_t0 * K_ygama0; % This is the camber aligning stiffness
    
    
end





