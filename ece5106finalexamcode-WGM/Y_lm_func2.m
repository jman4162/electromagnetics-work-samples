function [Y_lm] = Y_lm_func2(L, m, theta, phi)

x = cos(theta);
P_lm = legendre(L, x);
P_lm = P_lm(m+1);

Y_lm = sqrt((2.*L+1)./(4.*pi).*factorial(L-m)./factorial(L+m)).*P_lm.*exp(1i.*m.*phi);

end