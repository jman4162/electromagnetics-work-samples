% John Hodge
% ECE 5106 HW 5 P5c

clear all; clc; close all;

%% Define Constants
r = 0.6;
k = 1;
L = 1;
M = 1;
h = 1e-5;
theta = pi/4;
phi = pi/4;
a = 0.6;

x = k*r;

%% Define position
xx = 1;

x_arr = linspace(-1, 1, 50);
z_arr = linspace(0, 2, 50);

E_tot_mag_arr = zeros(length(x_arr), length(z_arr));

for xx = 1:length(x_arr)
    for zz = 1:length(z_arr)
        xx
        
        x = x_arr(xx);
        y = 0;
        z = z_arr(zz);
        
        r = sqrt(x.*x + y.*y + z.*z);
        theta = atan(sqrt(x.*x + y.*y) ./ z);
        phi = atan(y./x);
        
        theta_d = (180/pi) * theta;
        phi_d = (180/pi) * phi;
        
        
        %% Initialize Variables
        
        % Set field components to zero
        E_tot_phi = 0;
        E_tot_theta = 0;
        E_tot_r = 0;
        
        if (r > 0.6)
            
            L_term = 10;
            for L = 1:L_term
                
                % %% Define Y_LM
                % Y_LM = func_calc_Y_LM(L,M,theta,phi);
                
                %% Spherical Harmonics
                
                % Calc Y_LM and theta/phi derivatives
                Y_LM = func_calc_Y_LM(L,M,theta,phi);
                Y_LM_arr(1,xx) = Y_LM;
                
                dY_LM_dtheta = func_calc_dY_LM_dtheta(L,M,theta,phi,h);
                dY_LM_dtheta_arr(1,xx) = dY_LM_dtheta;
                
                dY_LM_dphi = func_calc_dY_LM_dphi(L,M,theta,phi,h);
                dY_LM_dphi_arr(1,xx) = dY_LM_dphi;
                
                % Define Hankel Function of First Kind
                h_L_1 = function_hankel_first(x,L);
                
                H = sqrt(pi/(2*x))*besselh(L+1/2,x);
                
                % Define Spherical Bessel Function
                j_L_kr = function_spherical_bessel(x,L);
                
                % Calculate alpha coeff
                xa = k*a;
                alpha_L = -(2*function_spherical_bessel(xa,L))./function_hankel_first(xa,L);
                
                dxjl_dx = ( (xa+h)*function_spherical_bessel((xa+h),L) - (xa-h)*function_spherical_bessel((xa-h),L) ) ./ (2*h);
                
                dxhl_dx = ( (xa+h)*function_hankel_first((xa+h),L) - (xa-h)*function_hankel_first((xa-h),L) ) ./ (2*h);
                
                beta_L = -2*dxjl_dx/dxhl_dx;
                
                % Calculate g and f functions
                x = k*r;
                g_L_kr = j_L_kr + (1/2) .* alpha_L .* h_L_1;
                
                f_L_kr = j_L_kr + (1/2) .* beta_L .* h_L_1;
                
                f_L_kr_test = function_calc_f_L_kr(x,L,h);
                
                drflr_dr = ( (x+h)*function_calc_f_L_kr((x+h),L,h) - (x-h)*function_calc_f_L_kr((x-h),L,h) ) ./ (2.*h);
                
                %% Break Up Equations into Terms
                
                % Term 1 components
                term1_phi = g_L_kr * (-1i./sqrt(L.*(L+1))) .* dY_LM_dtheta;
                
                term1_theta = g_L_kr * (-1i./sqrt(L.*(L+1))) .* (-1./sin(theta)) .* dY_LM_dphi;
                
                % Term 2 components
                term2_r = (1i).*sqrt(L.*(L+1))./(k*r) .* f_L_kr .* Y_LM;
                
                % Term 3 components
                term3_theta = 1./(k*r) .* drflr_dr .* (-1i./sqrt(L.*(L+1))) .* -dY_LM_dtheta;
                
                term3_phi = 1./(k*r) .* drflr_dr .* (-1i./sqrt(L.*(L+1))) .* -(1./sin(theta)).*dY_LM_dphi;
                
                %% Put Terms into E-field Components
                
                % E-Field Components L terms
                E_tot_phi_term = (1i).^L .* sqrt(4*pi*(2*L+1)) .* ( term1_phi + term3_phi );
                
                E_tot_theta_term = (1i).^L .* sqrt(4*pi*(2*L+1)) .* ( term1_theta + term3_theta );
                
                E_tot_r_term = (1i).^L .* sqrt(4*pi*(2*L+1)) .* ( term2_r );
                
                % Sum E-field L Terms
                E_tot_phi = E_tot_phi + E_tot_phi_term;
                E_tot_theta = E_tot_theta + E_tot_theta_term;
                E_tot_r = E_tot_r + E_tot_r_term;
                
            end
            
            E_tot_phi;
            E_tot_theta;
            E_tot_r;
            
            E_tot_mag_B = norm([E_tot_phi E_tot_theta E_tot_r]);
            
            E_tot_mag = sqrt( abs(E_tot_phi)^2 + abs(E_tot_theta)^2 + abs(E_tot_r)^2 );
            
        end
        
        if (r <= 0.6)
            E_tot_mag = 0;
        end
        
        E_tot_mag_arr(xx,zz) = E_tot_mag;
        
    end
end

figure;
imagesc(x_arr, z_arr, 20.*log10(E_tot_mag_arr'))
colorbar;
xlabel('x (m)')
ylabel('z (m)')
title('HW5 P5c: |E_{tot}| w/ k=1m^{-1}')
