% John Hodge
% ECE 5105 Final Exam Problem 2a
% 12/12/15
% TM Dielectric Slab Wave Guide

clear all; clc; close all;


A_cl_1 = 1; % A/m
d = 0.5 * 10^-6; %m
mu_0 = 4*pi*10^-7; % H/m
lambda = 1.0 * 10^-6; % m
c = 3*10^8; % m/s

w = (2*pi*(c/lambda));

eps_0 = 8.854 * 10^-12; % F/m
eps_co = 9*eps_0;
eps_cl_1 = 2*eps_0;
eps_cl_2 = 4*eps_0;

beta_i = sqrt(w^2*eps_cl_2*mu_0);
beta_f = sqrt(w^2*eps_co*mu_0);
beta_pts = 1000;

T_TM_11_arr = zeros(1, beta_pts);

x = 0;
for beta = linspace(beta_i, beta_f, beta_pts)
    x = x + 1;
    
    kappa_x1 = sqrt(beta^2 - w^2*eps_cl_1*mu_0);
    k_x_i = sqrt(w^2*eps_cl_1*mu_0 - beta^2);
    
    k_x_co = sqrt(w^2*eps_co*mu_0 - beta^2);
    
    kappa_x2 = sqrt(beta^2 - w^2*eps_cl_2*mu_0);
    k_x_f = sqrt(w^2*eps_cl_2*mu_0 - beta^2); 
    
    D_i_TM = [[1,1]; [k_x_i/eps_cl_1, -k_x_i/eps_cl_1]]; 
    D_t_TM = [[1,1]; [k_x_f/eps_cl_2, -k_x_f/eps_cl_2]];
    D_co_TM = [[1,1]; [k_x_co./eps_co, -k_x_co./eps_co]];
    
    D_i_TM_inv = inv(D_i_TM);
    D_t_TM_inv = inv(D_t_TM);
    D_co_TM_inv = inv(D_co_TM);
    
    P_co_TM = [[exp(-1i*k_x_co*d), 0]; [0, exp(+1i*k_x_co*d)]];
    
    T_TM = [D_t_TM_inv*D_co_TM*P_co_TM]*[D_co_TM_inv*D_i_TM];
    
    T_TM_11 = T_TM(1,1);
    T_TM_11_arr(x) = T_TM_11;
    
    T_TM_inv = inv(T_TM);
    
    
end

beta_plot = linspace(beta_i, beta_f, beta_pts);

figure; 
plot(beta_plot, abs(T_TM_11_arr))
xlabel('Beta (1/m)', 'Fontsize', 16)
ylabel('abs((T_{TM})_{11})', 'Fontsize', 16)
title('(T_{TM})_{11} vs. Progation Constant', 'Fontsize', 16)
grid on;
set(gca, 'FontSize', 16)

if(max(abs(T_TM_11_arr)) > 10)
   axis([beta_i beta_f 0 10])
end

figure; 
semilogy(beta_plot, abs(T_TM_11_arr))
xlabel('Beta (1/m)', 'Fontsize', 16)
ylabel('abs((T_{TM})_{11})', 'Fontsize', 16)
title('Log Scale: (T_{TM})_{11} vs. Progation Constant', 'Fontsize', 16)
grid on;
set(gca, 'FontSize', 16)