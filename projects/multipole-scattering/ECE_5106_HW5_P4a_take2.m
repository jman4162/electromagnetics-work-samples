
clear all; clc; close all;

%% Define Constants

k = 1;
% r = 1;
h = 1e-5;

plot_on = 0;

xSweep = 0:0.1:1;
zSweep = 0:0.1:1;

exEinc_re_arr = zeros(length(xSweep), length(zSweep));

for xx = 1:length(xSweep)
    for zz = 1:length(zSweep)
    xx   

x = xSweep(xx);
y = 0;
z = zSweep(zz);

r = sqrt(x.*x + y.*y + z.*z);

theta = atan(z./sqrt(x.*x+y.*y));

phi = atan(y./x);

% % Define direction
% phi_d = 0;
% theta_d = 90;
% 
% % Convert Angles from Degrees to Radians
% phi = (pi/180).*phi_d;
% theta = (pi/180).*theta_d;


exEinc = 0;
exEinc_mag_byTerm = zeros(1,20);
for Lx = 1:20
    
    % Define Bessel/Hankel Function Order
    L = Lx;
    m = 1;
    
    % Spherical Harmonic
    Y_lm = Y_lm_func2(L, m, theta, phi);
    
    % Spherical Harmonic First Derivative
    dY_lm_dtheta = (Y_lm_func2(L, m, theta + h, phi) - Y_lm_func2(L, m, theta - h, phi)) ./ (2.*h);
    
    dY_lm_dphi = (1i.*m).*Y_lm_func2(L, m, theta, phi);
    
    %% Calculate LY_phi
    
    exXlm_phi = -1i./sqrt(L.*(L+1)).*(-sin(theta).*dY_lm_dtheta);
    
    exXlm_theta = -1i./sqrt(L.*(L+1)).*(-cos(theta).*cos(phi).*(1./sin(theta)).*dY_lm_dphi);
    
    exXlm = -1i./sqrt(L.*(L+1)).*(-sin(theta) .* dY_lm_dtheta - cos(theta).*cos(phi) .* (1./sin(theta)) .* dY_lm_dphi);
    
    ex_d_er_cr_Xl1 = -1i./sqrt(L.*(L+1)).*(-cos(theta).*cos(phi).*dY_lm_dtheta + dY_lm_dphi);
    
    drjl_dr = ( ((r+h).*besselj_sph(L, k.*(r+h))) - ((r-h).*besselj_sph(L, k.*(r-h))) ) ./ (2.*h);
    
    r_plot = linspace(0.0001, 10, 1000);
    
    drjl_dr_plot = ( ((r_plot+h).*besselj_sph(L, k.*(r_plot+h))) - ((r_plot-h).*besselj_sph(L, k.*(r_plot-h))) ) ./ (2.*h);
    
    if (plot_on == 1)
    figure;
    hold on;
    plot(r_plot, r_plot.*besselj_sph(L, k.*r_plot))
    hold on;
    plot(r_plot, drjl_dr_plot)
    grid on;
    end
    
    ex_del_cr_jlXl1 = (1i.*sin(theta).*cos(phi).*sqrt(L.*(L+1))./r .* besselj_sph(L, k.*r).*Y_lm + (1./r).* drjl_dr .* ex_d_er_cr_Xl1);
    
    exEinc_term = (1i).^L.*sqrt(4.*pi.*(2.*L+1)) .* (besselj_sph(L, k.*r).* exXlm + (1/k).* ex_del_cr_jlXl1);
    exEinc_mag_byTerm(1, Lx) = abs(exEinc_term);
    
    exEinc = exEinc + exEinc_term;
    
end

figure(200);
hold on;
plot(1:20, exEinc_mag_byTerm); hold on;

exEinc_re_arr(xx,zz) = real(exEinc);

exEinc_re = real(exEinc);

    end
end

x
y
z
r
theta
phi

figure;
imagesc(zSweep, xSweep, exEinc_re_arr);