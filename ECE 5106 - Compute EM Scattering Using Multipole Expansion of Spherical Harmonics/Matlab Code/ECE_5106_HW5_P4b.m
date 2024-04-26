
clear all; clc; close all;

h = 1e-5;
k = 1; % k = 1m^-1;
r = 0.5;

theta_arr = linspace(0, pi, 100);
phi_arr = linspace(0, pi, 100);

ex_Einc = 0;

x_arr = linspace(0.0001, 6.*pi, 100);
y_arr = linspace(0.0001, 6.*pi, 100);

ex_Einc_re_arr = zeros(length(x_arr), length(y_arr));

for xx = 1:length(x_arr)
    for yy = 1:length(y_arr)
    xx
    
    x = x_arr(xx);
    y = y_arr(yy);
    z = 0;
    
    r = sqrt(x.*x + y.*y + z.*z);
    theta = atan(sqrt(x.*x + y.*y) ./ z);
    phi = atan(y./x);
    
    theta_d = (180/pi) * theta;
    phi_d = (180/pi) * phi;
    
    ex_Einc = 0;
    L_num = 35;
    for Lx = 1:L_num
        
        % Define Parameters
        L = Lx;
        M = 1;
        %phi = phi_arr(x);
        %theta = pi/2;
        
        % ex dot er, etheta, and ephi
        exer = sin(theta).*cos(phi);
        exet = cos(theta).*cos(phi);
        exep = -sin(phi);
        
        % Define Spherical Bessel Function
        j_L_kr = function_spherical_bessel(k*r,L);
        
        % Calc d[r*j_l_kr]/dr
        drjl_dr = ( (r+h).*function_spherical_bessel(k*(r+h),L) - (r-h).*function_spherical_bessel(k*(r-h),L) ) ./ (2.*h);
        
        % Calc Y_LM and theta/phi derivatives
        Y_LM = func_calc_Y_LM(L,M,theta,phi);
        Y_LM_arr(1,xx) = Y_LM;
        
        dY_LM_dtheta = func_calc_dY_LM_dtheta(L,M,theta,phi,h);
        dY_LM_dtheta_arr(1,xx) = dY_LM_dtheta;
        
        dY_LM_dphi = func_calc_dY_LM_dphi(L,M,theta,phi,h);
        dY_LM_dphi_arr(1,xx) = dY_LM_dphi;
        
        % ex dot X_LM
        ex_d_Xlm = -1i./sqrt(L.*(L+1)) .* ( -sin(phi)*dY_LM_dtheta - cos(theta)*cos(phi)*(1./sin(theta)).*dY_LM_dphi );
        
        % ex dot (er cross X_LM)
        ex_d_er_cr_Xlm = -1i./sqrt(L*(L+1)) .* ( -cos(theta)*cos(phi)*dY_LM_dtheta + dY_LM_dphi );
        
        term1 =  j_L_kr .* ex_d_Xlm;
        
        term2 = 1i.*sqrt(L*(L+1))./r .* j_L_kr .* Y_LM .* sin(theta).*cos(phi);
        
        term3 = (1/r) * drjl_dr * ex_d_er_cr_Xlm;
        
        ex_Einc_term = (1i).^L .* sqrt(4*pi*(2*L+1)) .* ( term1 + (1/k)*(term2 + term3) );
        
        ex_Einc_term_L_mag(1,Lx) = abs(ex_Einc_term);
        
        ex_Einc = ex_Einc + ex_Einc_term;
        
        
    end
    ex_Einc_re = real(ex_Einc);
    ex_Einc_re_arr(xx, yy) = ex_Einc_re;
    % Show magnitude of each L term
    figure(200);
    hold on;
    plot(1:L_num, ex_Einc_term_L_mag); hold on;
    
    end
end

% figure;
% hold on;
% plot(phi_arr, Y_LM_arr); hold on;
% %plot(theta_arr, dY_LM_dtheta_arr); hold on;
% plot(phi_arr, dY_LM_dphi_arr); hold on;
% grid on;

figure;
imagesc(x_arr, y_arr, ex_Einc_re_arr')
xlabel('x (m)')
ylabel('y (m)')
title('Re[ex dot E_{inc}] w/ k=1m^{-1}')
colorbar;