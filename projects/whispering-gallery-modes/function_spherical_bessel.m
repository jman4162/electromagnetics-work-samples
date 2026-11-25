% defines spherical bessel function

function fval = function_spherical_bessel(x,L)

fval = sqrt(pi/2./x) .* besselj(L+1/2,x);