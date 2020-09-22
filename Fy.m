function [ Fy ] = Fy(alpha,Fz,gamma)

Cy = p_cy1*lambda_Cy;
Muy = (p_Dy1 + p_Dy2*df_zfr)*(1-(p_Dy3*(sin(gamma))^2)*lambda_muV);
Dy = Muy*Fz;
Kyalpha0 = p_Ky1*Fz*lambda_Fz0*sin(2*atan(Fz/(p_Ky2*Fz*lambda_Fz0)))*lambda_Ey;
Kyalpha = Kyalpha0*(1-p_Ky3*(sin(gamma))^2);
By = Kyalpha/(Cy*Dy+e);
SHy = (p_Hy1 + p_Hy2*df_zfr)*lambda_muxy+p_Hy3*sin(gamma)*lambda_Kygamma;
alphay = alpha + SHy;
SVy = DF_zfr*((p_Vy1+p_Vy2*df_zfr)*lambda_muxy+(p_Vy3+p_Vy4*df_zfr)*sin(gamma)*lambda_Kygamma)*lambda_muV;
Ey = (p_Ey1 + p_Ey2*df_zfr)*(1-(p_Ey3 + p_Ey4*sin(gamma))*sign(alphay))*lambda_Ey;
Kygamma0 = (p_Hy3*K_yalpha0+Fz*(p_V3_p_Vy3*df_zfr))*lambda_Kygamma;

Fy = Dy*sin(Cy*atan(By*alphay-Ey(By*alphay-atan(By*alphay))))+SVy;

end