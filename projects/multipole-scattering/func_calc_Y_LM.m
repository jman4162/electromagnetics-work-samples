% This function evaluate spherical harmonics function
% Y_lm(theta,phi)
% refer to the solution

function val = func_calc_Y_LM(L,M,theta,phi)

x = cos(theta);

factor_1_array = legendre(L,x,'norm');

factor_1 = factor_1_array(M+1);

factor_2 = (-1)^M * exp(1i*M*phi) / sqrt(2*pi);

val = factor_1 .* factor_2;