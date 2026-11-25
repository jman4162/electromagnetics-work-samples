function f_L_kr = function_calc_f_L_kr(x,L,h)

%% Define Hankel Function of First Kind
h_L_1 = function_hankel_first(x,L);

H = sqrt(pi/(2*x))*besselh(L+1/2,x);

%% Define Spherical Bessel Function
j_L_kr = function_spherical_bessel(x,L);

%% Calculate alpha coeff
% alpha_L = -(2*j_L_kr)./h_L_1;

dxjl_dx = ( (x+h)*function_spherical_bessel((x+h),L) - (x-h)*function_spherical_bessel((x-h),L) ) ./ (2*h);

dxhl_dx = ( (x+h)*function_hankel_first((x+h),L) - (x-h)*function_hankel_first((x-h),L) ) ./ (2*h);

beta_L = -2*dxjl_dx/dxhl_dx;

%% Calculate g and f functions

%g_L_kr = j_L_kr + (1/2) .* alpha_L .* h_L_1

f_L_kr = j_L_kr + (1/2) .* beta_L .* h_L_1;

end