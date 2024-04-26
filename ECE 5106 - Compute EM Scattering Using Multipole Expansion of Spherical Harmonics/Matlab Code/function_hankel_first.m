function [h_L_1] = function_hankel_first(x,L)

%% Define Hankel Function of First Kind
%h_L_1 = (-1i).^(L+1) .* exp(1i.*x) ./ x;

h_L_1 = sqrt(pi/(2.*x)).*(besselj(L+1/2,x) + 1i.*bessely(L+1/2,x));

end