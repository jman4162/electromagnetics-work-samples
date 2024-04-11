
clear all; close all; clc;

phi_plot = 0.1:0.4:2*pi;
theta_plot = linspace(0.1,pi,250);

delta = 1e-7;

XLM_theta_arr = zeros(length(phi_plot),length(theta_plot));
XLM_phi_arr = zeros(length(phi_plot),length(theta_plot));
XLM_mag_arr = zeros(length(phi_plot),length(theta_plot));

for aa = 1:length(phi_plot);
    for bb = 1:length(theta_plot);
        
        phi = phi_plot(aa);
        theta = theta_plot(bb);
        
        L = 110;
        
        M = L;
        XLM_theta = function_XLM_theta_mod(theta,phi,L,M);
        XLM_phi = function_XLM_phi_mod(theta,phi,L,M,delta);
        
        XLM_theta_arr(aa,bb) = abs(XLM_theta/10e100);
        XLM_phi_arr(aa,bb) = abs(XLM_phi/10e100);
        XLM_mag_arr(aa,bb) = sqrt(abs(XLM_theta/10e100)^2 + abs(XLM_phi/10e100)^2);
        
    end
end

figure;
subplot(1,2,1)
imagesc(phi_plot, theta_plot, XLM_theta_arr')
colorbar;
xlabel('Phi');
ylabel('Theta');
title('X_{LM theta}')

subplot(1,2,2)
imagesc(phi_plot, theta_plot, XLM_phi_arr')
colorbar;
xlabel('Phi');
ylabel('Theta');
title('X_{LM phi}')

figure;
imagesc(phi_plot, theta_plot, XLM_mag_arr')
colorbar;
xlabel('Phi');
ylabel('Theta');
title('X_{LM mag}')

figure;
plot(theta_plot, XLM_mag_arr(1,:)')
xlabel('Theta');
ylabel('X_{LM mag}');
title('X_{LM mag} vs. Theta')

% figure;
% subplot(1,2,1)
% imagesc(XLM_theta_arr')
% colorbar;
% xlabel('Phi');
% ylabel('Theta');
% title('X_{LM theta}')
% 
% subplot(1,2,2)
% imagesc(XLM_phi_arr')
% colorbar;
% xlabel('Phi');
% ylabel('Theta');
% title('X_{LM phi}')


