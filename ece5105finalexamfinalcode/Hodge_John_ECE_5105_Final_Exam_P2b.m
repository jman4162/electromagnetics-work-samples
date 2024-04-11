% John Hodge
% ECE 5105 Problem 2b

clear all; clc; close all;

A_cl_1 = 1; % A/m
d = 0.5 * 10^-6; %m
mu_0 = 4*pi*10^-7; % H/m
lambda = 1.0 * 10^-6; % m
c = 3*10^8; % m/s

w = (2*pi*(c/lambda))

eps_0 = 8.854 * 10^-12; % F/m
eps_co = 9*eps_0;
eps_cl_1 = 2*eps_0;
eps_cl_2 = 4*eps_0;
z = 0;
t = 0;
x_1 = 0;
x_2 = d;

beta = 1.7971 * 10^7 % 1/m

kappa_x1 = sqrt(beta^2 - w^2*eps_cl_1*mu_0)
k_x_i = sqrt(w^2*eps_cl_1*mu_0 - beta^2);

k_x_co = sqrt(w^2*eps_co*mu_0 - beta^2)

kappa_x2 = sqrt(beta^2 - w^2*eps_cl_2*mu_0)
k_x_f = sqrt(w^2*eps_cl_2*mu_0 - beta^2);

A_co_1 = 1/2 + 1i.*(9/4).*(kappa_x1/k_x_co)

A_co_2 = 1/2 - 1i.*(9/4).*(kappa_x1/k_x_co)

A_cl_2 = A_co_1.*exp(-1i.*k_x_co.*d) + A_co_2.*exp(+1i.*k_x_co.*d)

x_pts = 30000;

x = linspace(-0.5, 1, x_pts);

x = x.*(10^-6);

H_cl_1 =  A_cl_1 .* exp(1i.*(w*t - beta*z)) .* exp(kappa_x1.*(x-x_1));

H_co = (A_co_1 .* exp(-1i.* k_x_co .* x) + A_co_2 .* exp(+1i.* k_x_co .* x)) .* exp(1i.*(w*t - beta*z));

H_cl_2 =  A_cl_2 .* exp(1i.*(w*t - beta*z)) .* exp(-kappa_x2.*(x-x_2));

H_plot = [H_cl_1(1:x_pts./3) H_co(x_pts./3 + 1:2.*x_pts./3) H_cl_2(2.*x_pts./3 + 1:end)];

figure;
hold on;
plot(x, H_plot)
grid on;
xlabel('Position on x-axis (m)', 'FontSize', 16)
ylabel('H-Field Strength (A/m)', 'FontSize', 16)
title('H_{y} Field Strength in Dielectric Slab Waveguide vs. Position (z=0; t=0)', 'FontSize', 16)
set(gca, 'FontSize', 16)
hold off;

hold on;
plot([x_1 x_1], [0 14]);
hold off;

hold on;
plot([x_2 x_2], [0 14]);
hold off;


S_z_1 = (beta./(2.*w.*eps_cl_1)).*real(A_cl_1).^2 .* exp(+2.*kappa_x1.*(x-x_1));



A_co_comb = (A_co_1 .* exp(-1i.* k_x_co .* x) + A_co_2 .* exp(+1i.* k_x_co .* x));
A_co_comb_conj = conj(A_co_comb);

S_z_2 = (beta./(2.*w.*eps_cl_1)) .* real(A_co_comb .* A_co_comb_conj);

S_z_3 = (beta./(2.*w.*eps_cl_2)).*real(A_cl_2).^2 .* exp(-2.*kappa_x1.*(x-x_2));

S_z = [S_z_1(1:x_pts./3) S_z_2(x_pts./3 + 1:2.*x_pts./3) S_z_3(2.*x_pts./3 + 1:end)];

% figure;
% plot(x(1:x_pts./3), S_z_1(1:x_pts./3))
% 
% figure;
% plot(x(2.*x_pts./3 + 1:end), S_z_3(2.*x_pts./3 + 1:end))

figure;
hold on;
plot(x, S_z)
grid on;
xlabel('Position on x-axis (m)', 'FontSize', 16)
ylabel('Poynting Vector Strength (W/m)', 'FontSize', 16)
title('Z-Component of Time-Averaged Poynting Vector vs. Position', 'FontSize', 16)
set(gca, 'FontSize', 16)
hold off;

hold on;
plot([x_1 x_1], [0 max(S_z)]);
hold off;

hold on;
plot([x_2 x_2], [0 max(S_z)]);
hold off;
