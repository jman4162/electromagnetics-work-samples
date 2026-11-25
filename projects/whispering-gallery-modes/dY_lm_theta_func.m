function dY_lm_theta = dY_lm_theta_func(L,m,theta_d,phi_d, h)

dY_lm_theta = (Y_lm_func(L,m,theta_d+h,phi_d) - Y_lm_func(L,m,theta_d-h,phi_d))./(2.*h);

end