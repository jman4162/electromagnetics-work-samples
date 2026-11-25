
clear all; clc; close all;

a = 0.6;
k = 1;
xa = k*a;
h = 1e-5;

L_alpha_beta_coeff_arr = zeros(6,3);

for L = 1:10;
    
    % Define Spherical Bessel Function
    j_L_kr = function_spherical_bessel(xa,L)
    
    % Define Spherical Hankel Function
    h_L_1 = function_hankel_first(xa,L)
    
    % Calc Alpha and Beta
    alpha_L = -2*j_L_kr./h_L_1
    
    dxjl_dr = ( (xa+h)*function_spherical_bessel(xa+h,L) - (xa-h)*function_spherical_bessel(xa-h,L) ) ./ (2*h)
    
    dxhl_dr = ( (xa+h)*function_hankel_first(xa+h,L) - (xa-h)*function_hankel_first(xa-h,L) ) ./ (2*h)
    
    beta_L = -2*dxjl_dr./dxhl_dr
    
    g_L_kr = j_L_kr + (1/2) .* alpha_L .* h_L_1
    
    f_L_kr = j_L_kr + (1/2) .* beta_L .* h_L_1
    
    f_L_kr_test = function_calc_f_L_kr(xa,L,h)
    
    f_L_kr_test2 = function_calc_f_L_kr_B(xa,L,h)
    
    dxf_L_kr_dr = ( (xa+h)*function_calc_f_L_kr_B(xa+h,L,h) - (xa-h)*function_calc_f_L_kr_B(xa-h,L,h) ) ./ (2*h)
    
    L_alpha_beta_coeff_arr(L,1) = L;
    L_alpha_beta_coeff_arr(L,2) = alpha_L;
    L_alpha_beta_coeff_arr(L,3) = beta_L;
end