% calculate the theta component of a vector spherical harmonic

function fval = function_XLM_theta_mod(theta,phi,L,M)

x = cos(theta);

norm_const = sqrt((2*L+1) / 4 / pi );

PLM_array = legendre(L,x);

PLM = PLM_array(M+1);

fval = - M / sin(theta) * norm_const * PLM * exp(1j*M*phi) / sqrt(L*(L+1));


