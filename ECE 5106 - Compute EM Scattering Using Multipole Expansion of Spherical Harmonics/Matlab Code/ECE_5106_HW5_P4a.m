
clear all; clc; clc;

k = 1;

phi_d = 0;
theta_d = 90;

phi = phi_d.*pi./180;
theta = theta_d.*pi./180;


L = 1;
m = 1;

r = 1;
G = besselj_sph(L, k.*r);

[Y_lm] = Y_lm_func2(L, m, theta, phi)

h = 1e-4;

% Convert Angles from Degrees to Radians
phi = (pi/180).*phi_d;
theta = (pi/180).*theta_d;

% Spherical Harmonic
Y_lm = Y_lm_func2(L, m, theta, phi);

% Spherical Harmonic First Derivative
dY_lm_dtheta = (Y_lm_func2(L, m, theta + h, phi) - Y_lm_func2(L, m, theta - h, phi)) ./ (2.*h);

%% Vector Spherical Harmonics X_lm(theta,phi)


