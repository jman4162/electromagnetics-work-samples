
clear all; clc; close all;

a = 0.6;
k = 1;
xa = k*a;
h = 1e-5;

load L_alpha_beta_coeff_arr;

term_sum = 0;
L_terms = 6;
term_mag_arr = zeros(1, L_terms);
for L = 1:L_terms
   
    L;
    alpha_L = L_alpha_beta_coeff_arr(L,2);
    beta_L = L_alpha_beta_coeff_arr(L,3);
    
    term = (2*L+1) .* ( abs(alpha_L)^2 + abs(beta_L)^2 );
    
    term_mag_arr(1,L) = abs(term);
    
    term_sum = term_sum + term;
    
    
end

constant = pi./(2*k*k);

rcs_sc = constant * term_sum

figure;
subplot(2,1,1)
plot(1:L_terms, term_mag_arr)
xlabel('L')
ylabel('Scattering L-term Magnitude')
title('RCS Scattering L-term Magnitude vs. L')
grid on;
subplot(2,1,2)
plot(1:L_terms, 10.*log10(term_mag_arr))
title('dB')
grid on;
xlabel('L')
ylabel('Scattering L-term Magnitude')