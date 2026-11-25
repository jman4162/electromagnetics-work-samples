function [Y_lm] = Y_lm_func(L,m,theta_d,phi_d)

%% Calculate Spherical Harmonic
phi = phi_d.*(pi./180);
theta = theta_d.*(pi./180);
x = cosd(theta_d);
P_l_m = legendre(L,x);
P_l_10 = P_l_m(m+1);
Y_l_m = sqrt((2.*L+1)./(4.*pi).*factorial(L-m)./factorial(L+m)).*P_l_10.*exp(1i.*m.*phi);

Y_lm = Y_l_m;

end