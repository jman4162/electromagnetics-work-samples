% John Hodge
% ECE 5105 Final Exam Problem 1
% 12/12/15
% TE Three layer film

clear all; clc; close all;

for z = 1:2
    
    d_1 = 0.01;
    d_step = 0.01;
    d_fin = 10;
    
    d_len = length(d_1:d_step:d_fin);
    
    R_arr = zeros(1, d_len);
    T_arr = zeros(1, d_len);
    
    x = 0;
    for d = d_1:d_step:d_fin; % distance in mm
        x = x + 1;
        d = 10^-3*d;
        
        theta_i = 60;
        eps_0 = 8.854*10^-12;
        mu_0 = 4*pi*10^-7;
        
        w = 3*10^11;
        
        eps_1 = 4.*eps_0;
        eps_3 = 4.*eps_0;
        
        eps_2A = 1.0*eps_0;
        eps_2B = 8.0.*eps_0;
        eps_2X = [eps_2A, eps_2B];
        
        eps_2 = eps_2X(z);
        
        k = w.*sqrt(eps_1.*mu_0);
        k_1x = (k)./sqrt(1+tand(theta_i).^2);
        
        beta = tand(theta_i).*k_1x;
        
        k_2x = sqrt(w.^2.*eps_2.*mu_0 - beta.^2);
        k_3x = sqrt(w.^2.*eps_3.*mu_0 - beta.^2);
        
        D_1_TE = [[1,1];[k_1x/mu_0, -k_1x/mu_0]];
        D_2_TE = [[1,1];[k_2x/mu_0, -k_2x/mu_0]];
        D_3_TE = [[1,1];[k_3x/mu_0, -k_3x/mu_0]];
        
        D_2_TE_inv = inv(D_2_TE);
        D_3_TE_inv = inv(D_3_TE);
        
        P_2_TE = [[exp(-1i.*k_2x.*d), 0]; [0, exp(+1i.*k_2x.*d)]];
        
        T_TE = [D_3_TE_inv*D_2_TE*P_2_TE]*[D_2_TE_inv*D_1_TE];
        T_TE_inv = inv(T_TE);
        T_TE_inv_11 = T_TE_inv(1,1);
        T_TE_inv_21 = T_TE_inv(2,1);
        
        
        R = abs(T_TE_inv_21/T_TE_inv_11).^2;
        R_arr(1,x) = R;
        
        T = (mu_0*k_3x)./(mu_0*k_1x).*1./abs(T_TE_inv_11).^2;
        T_arr(1,x) = T;
        
        
    end
    
    
    d_plot = d_1:d_step:d_fin;
    
    if (z == 1)
        a = 1;
        b = 3;
    elseif (z == 2)
        a = 2;
        b = 4;
    end
    
    figure(100);
    subplot(2,2,a)
    plot(d_plot, R_arr)
    axis([d_1 d_fin 0 1])
    title('Reflected Power vs. Layer Thickness', 'Fontsize', 16)
    xlabel('Distance (mm)', 'Fontsize', 16)
    ylabel('Reflected Power (R)', 'Fontsize', 16)
    set(gca, 'FontSize', 16)
    
    subplot(2,2,b)
    plot(d_plot, T_arr)
    axis([d_1 d_fin 0 1])
    title('Transmitted Power vs. Layer Thickness', 'Fontsize', 16)
    xlabel('Distance (mm)', 'Fontsize', 16)
    ylabel('Transmitted Power (T)', 'Fontsize', 16)
    set(gca, 'FontSize', 16)
    
    figure(200)
    subplot(2,1,a)
    plot(d_plot, R_arr, 'r--', d_plot, T_arr, 'b-')
    axis([d_1 d_fin 0 1])
    if (z == 1)
        title('Case A: Reflected Power vs. Layer Thickness', 'Fontsize', 16)
    elseif (z == 2)
        title('Case B: Reflected Power vs. Layer Thickness', 'Fontsize', 16)
    end
    xlabel('Distance (mm)', 'Fontsize', 16)
    ylabel('Reflected Power (R)', 'Fontsize', 16)
    set(gca, 'FontSize', 16)
    legend('Reflected', 'Transmitted')
    grid on;
    
    figure(300)
    subplot(2,1,a)
    plot(d_plot, 10*log10(R_arr), 'r--', d_plot, 10*log10(T_arr), 'b-')
    %axis([d_1 d_fin 0 1])
    if (z == 1)
        title('Case A: Reflected Power vs. Layer Thickness', 'Fontsize', 16)
    elseif (z == 2)
        title('Case B: Reflected Power vs. Layer Thickness', 'Fontsize', 16)
    end
    xlabel('Distance (mm)', 'Fontsize', 16)
    ylabel('Reflected Power (R)', 'Fontsize', 16)
    set(gca, 'FontSize', 16)
    legend('Reflected', 'Transmitted')
    grid on;
    
end