
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
w = 1.970343823438235e+15; % Set Resonant Frequency
L = 3;

r_plot = linspace(0, a, 1000);
k_co = (w/c)*n_co;
k_co_a = k_co*a;
Y = besselj_sph(L, k_co.*r_plot);

figure;
plot(r_plot, Y);