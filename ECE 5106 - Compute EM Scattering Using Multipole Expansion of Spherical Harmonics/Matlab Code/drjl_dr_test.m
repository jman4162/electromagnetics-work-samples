close all; clc; clear all;


% test this function 
    r_plot = linspace(0, 1, 100);
    drjl_dr_test = ( (r_plot+h).*function_spherical_bessel(k*(r_plot+h),L) - (r_plot-h).*function_spherical_bessel(k*(r_plot-h),L) ) ./ (2.*h);
    
    figure;
    plot(r_plot, drjl_dr_test)