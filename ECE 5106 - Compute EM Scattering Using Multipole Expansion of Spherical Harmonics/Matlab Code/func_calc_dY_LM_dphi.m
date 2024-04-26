function dY_LM_dphi = func_calc_dY_LM_dphi(L,M,theta,phi,h)

dY_LM_dphi = (func_calc_Y_LM(L,M,theta,phi+h) - func_calc_Y_LM(L,M,theta,phi-h))/(2.*h);

end


