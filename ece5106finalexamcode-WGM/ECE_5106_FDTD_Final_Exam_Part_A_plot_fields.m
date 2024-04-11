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

% Calc Freq Range
f_min = 3.059e14;
w_min = 2*pi*f_min;

f_max = 3.224e14;
w_max = 2*pi*f_max;

%w = 1.970343823438235e+15; % Set Resonant Frequency
%w = 1.924369773697737e+15;
% w = 1.928929639296393e+15; % L = 50;
% w = 1.934413844138442e+15; % L = 110;
%w = 1.970092290922909e+15; % L + 85;
%w = 1.966281062810628e+15; %L = 130;
%w = 2.009706473064731e+15; %L = 120;
%w = 1.997465574655747e+15; % L = 125;
w = 1.9652e15;

k_co = (w/c)*n_co;
k_cl = (w/c)*n_cl;

L = 130;

phi_plot = linspace(0.01,2*pi,3);
r_plot = linspace(0.01*a, a, 1000);

E_mag_arr = zeros(length(r_plot), length(phi_plot));

for r_x = 1:length(r_plot)
    for phi_x = 1:length(phi_plot)
        
        r = r_plot(r_x);
        phi = phi_plot(phi_x);
        theta = pi/2;
        %x = 2;
        
      %  Y1 = function_hankel_first(x,L);
        
        %% Define Spherical Bessel Function
        %Y2 = besselj_sph(L, x);
        M = L;
        XLM_theta = function_XLM_theta_mod(theta,phi,L,M)/10e120;
        XLM_phi = function_XLM_theta_mod(theta,phi,L,M)/10e120;
        
        Ch_Eq_L = (1./Z_co).*( besselj_sph(L-1, k_co.*a) - L/(k_co.*a).*besselj_sph(L, k_co.*a) ) ./ besselj_sph(L, k_co.*a);
        
        Ch_Eq_R = (1./Z_cl).*( function_hankel_first(L-1, k_cl.*a) - L/(k_cl.*a).*function_hankel_first(L, k_cl.*a) ) ./ function_hankel_first(L, k_cl.*a);
        
        y_test = Ch_Eq_L - Ch_Eq_R;
        
        %r = 0.8*a;
        
        E_theta = besselj_sph(L, k_co.*r)*XLM_theta;
        
        E_phi = besselj_sph(L, k_co.*r)*XLM_phi;
        
        E_mag = sqrt( abs(E_theta).^2 + abs(E_phi).^2 );
        
        E_mag_arr(r_x, phi_x) = E_mag;
        
    end
end

figure;
subplot(1,2,1)
imagesc(r_plot, phi_plot, E_mag_arr'/max(max(E_mag_arr)))
colorbar;
xlabel('r');
ylabel('Phi');
title('E_{mag} L = 130, M = 130, theta = pi/2 (normalized)')

subplot(1,2,2)
plot(r_plot, E_mag_arr(:,1)/max(E_mag_arr(:,1)));
xlabel('r (meters)')
ylabel('E-Field Magnitude (normalized)')
title('L = 130, M = 130, phi = 0, theta = pi/2')
grid on;