% calculate the phi component of a vector spherical harmonic

function fval = function_XLM_phi_mod(theta,phi,L,M,delta)

x = cos(theta);

x_p = cos(theta + delta);
x_m = cos(theta - delta);

norm_const = -1i * ...
    sqrt((2*L+1) / 4 / pi ) * ...
    exp(1i*M*phi);

PLM_array_m = legendre(L,x_m);  PLM_m = PLM_array_m(M+1);
PLM_array_p = legendre(L,x_p);  PLM_p = PLM_array_p(M+1);

dPLM_dtheta = (PLM_p - PLM_m) / 2 / delta;

fval = norm_const * dPLM_dtheta / sqrt(L*(L+1));