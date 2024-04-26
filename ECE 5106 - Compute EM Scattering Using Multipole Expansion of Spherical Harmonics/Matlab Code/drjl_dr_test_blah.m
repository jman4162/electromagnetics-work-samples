close all; clc; clear all;

h = 1e-5;
k = 1; % k = 1m^-1;
r = 0.5;

L = 1;

% test this function 
    r_plot = linspace(0, 10, 1000);
    drjl_dr_test = ( (r_plot+h).*function_spherical_bessel(k*(r_plot+h),L) - (r_plot-h).*function_spherical_bessel(k*(r_plot-h),L) ) ./ (2.*h);
    
    figure;
    hold on;
    plot(r_plot, (r_plot).*function_spherical_bessel(k*(r_plot),L)); hold on;
    plot(r_plot, drjl_dr_test); hold on;
    grid on;