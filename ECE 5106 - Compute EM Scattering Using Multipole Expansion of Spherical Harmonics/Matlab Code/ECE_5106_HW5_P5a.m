% ECE 5106 HW5 P5

clear all; clc; close all;

k = 1;
r = 0.6;
h = 1e-5;

x = k*r;

L_alpha_beta_coeff_arr = zeros(3,6);

for Lx = 1:1

L = Lx;

% Define Hankel Function of First Kind
% h_L_1 = (-1i).^(L+1) .* exp(1i.*x) ./ x;
h_L_1 = function_hankel_first(x,L)

% Define Spherical Bessel Function
j_L_kr = function_spherical_bessel(x,L);

% Calculate alpha coeff
alpha_L = -(2*j_L_kr)./h_L_1;

dxjl_dx = ( (x+h)*function_spherical_bessel((x+h),L) - (x-h)*function_spherical_bessel((x-h),L) ) ./ (2*h);

dxhl_dx = ( (x+h)*function_hankel_first((x+h),L) - (x-h)*function_hankel_first((x-h),L) ) ./ (2*h);

beta_L = -2*dxjl_dx/dxhl_dx

L_alpha_beta_coeff_arr(1,Lx) = Lx;
L_alpha_beta_coeff_arr(2,Lx) = alpha_L;
L_alpha_beta_coeff_arr(3,Lx) = beta_L;

end

L_alpha_beta_coeff_arr = L_alpha_beta_coeff_arr';
