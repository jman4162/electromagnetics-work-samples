function dY_LM_dtheta = func_calc_dY_LM_dtheta(L,M,theta,phi,h)

dY_LM_dtheta = (func_calc_Y_LM(L,M,theta+h,phi) - func_calc_Y_LM(L,M,theta-h,phi))/(2.*h);

end


