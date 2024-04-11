% John Hodge
% 05/03/16
% ECE 5106 Final Exam
% Spherical WGM Resonator

clear all; clc; close all;

%% Define constants and key values
n_co = 1.4;
n_cl = 1.1;

eps_r_co = sqrt(n_co);
eps_r_cl = sqrt(n_cl);

eps_0 = 8.854e-12;
mu_0 = 4*pi*1e-7;
c = 1/sqrt(eps_0*mu_0);

eps_co = eps_r_co*eps_0;
eps_cl = eps_r_cl*eps_0;

Z_0 = sqrt(mu_0/eps_0);
Z_co = sqrt(mu_0/eps_co);
Z_cl = sqrt(mu_0/eps_cl);

a = 15e-6;

y_test_min = 1e100;
L_y_min = 0;
w_y_min = 0;

% Calc Freq Range
f_min = 3.059e14;
w_min = 2*pi*f_min;

f_max = 3.224e14;
w_max = 2*pi*f_max;

for L = 30;
    phi = 0.001;
    theta = 0.001;
    x = 2;
    
    Y1 = function_hankel_first(x,L);
    
    %% Define Spherical Bessel Function
    Y2 = besselj_sph(L, x);
    M = L;
    XLM_theta = function_XLM_theta(theta,phi,L,M);
    XLM_phi = function_XLM_theta(theta,phi,L,M);
    
    
    n_pts = 100000;
    
    Ch_Eq_L_arr = zeros(1,n_pts);
    Ch_Eq_R_arr = zeros(1,n_pts);
    y_test_arr = zeros(1,n_pts);
    
    %w_plot = linspace(w_min, w_max, n_pts);
    w_plot = linspace(1.924e15, 1.925e15, n_pts);
    
    ii = 0;
    for ww = 1%:n_pts
        ii = ii + 1;
        
        %w = w_plot(ww);
        %w = 1.970343823438235e+15;
        w = 1.924369773697737e+15;
        
        lambda_min = 0.93e-6;
        lambda_max = 0.98e-6;
        
        k_co = (w/c)*n_co;
        k_cl = (w/c)*n_cl;
        
        Ch_Eq_L = (1./Z_co).*( besselj_sph(L-1, k_co.*a) - L/(k_co.*a).*besselj_sph(L, k_co.*a) ) ./ besselj_sph(L, k_co.*a);
        
        Ch_Eq_R = (1./Z_cl).*( function_hankel_first(L-1, k_cl.*a) - L/(k_cl.*a).*function_hankel_first(L, k_cl.*a) ) ./ function_hankel_first(L, k_cl.*a);
        
        y_test = Ch_Eq_L - Ch_Eq_R;
        
        Ch_Eq_L_arr(ii) = Ch_Eq_L;
        Ch_Eq_R_arr(ii) = Ch_Eq_R;
        y_test_arr(ii)  = y_test;
        
        if ((abs(y_test) < abs(y_test_min)) && (1.94e15 > w) ) % && (2.00e15 > w) )
            y_test_min = y_test;
            L_y_min = L;
            w_y_min = w;
        end
        
    end
    
    figure(200);
    hold on;
    plot(w_plot, Ch_Eq_L_arr); hold on;
    plot(w_plot, real(Ch_Eq_R_arr)); hold on;
    plot(w_plot, imag(Ch_Eq_R_arr)); hold on;
    ylim([-1 1])
    
    figure(300);
    hold on;
    plot(w_plot, y_test_arr); hold on;
    grid on;
    ylim([-1 1])
    
end
