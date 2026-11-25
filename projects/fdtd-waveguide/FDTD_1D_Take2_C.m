% John Hodge
% 05/01/16
% ECE 5106 1-D FDTD Code

clear all; clc; close all;

%% Set Constants
eps = 1;
mu = 1;

c = 1/sqrt(eps*mu);

L = 30;
%dz = 1;
N = 31;
dz = L./(N-1);
dt = 0.5*dz;

Mt  = 10000;

% t_plot = 

%% Step 1: Set Initial Values for Ex and Hy
Ex = zeros(N, Mt);
Hy = zeros(N, Mt);

% Ex(:,2) = sin(linspace(0, pi, N));
% Ex(:,2) = ones(1,N);

sig_terms = 5;
sig = ones(1, N-1);
for cc = 1:sig_terms
sig = sig + (1/cc).*sin((cc.*pi)./30.*(1:N-1));
end

figure;
plot(sig)

Ex(1:end,2) = sin((5*pi)./(N-1).*(1:N));
%Ex(1:end-1,2) = sig;
%Hy(1:end-1,2) = sig;
%Ex(:,2) = ones(1,N);
%Hy(:,2) = ones(1,N);
%Ex(:,2) = rand(1,N);

figure;
stem((1:N)*dz, Ex(:,2));
xlabel('Position on z-axis (norm units)', 'FontSize', 16)
ylabel('E_x Field Amplitude (V/m)', 'FontSize', 16)
title('1-D FDTD Initial Field Distribution', 'FontSize', 16)
grid on;
set(gca,'Fontsize',16);
xlim([0 30])

for m = 2:1:Mt
    
    % Update BCs
    Ex(1,m) = 0;
    Ex(N,m) = 0;
    Hy(N,m) = 0;
    
    % Update E 
    for n = 1:N-1
        Hy(n,m) = Hy(n,m-1) + dt./mu .* ( -(Ex(n+1,m) - Ex(n,m))./dz );
    end
    
        % Update BCs
    Ex(1,m) = 0;
    Ex(N,m) = 0;
    Hy(N,m) = 0;
    
    % Update H
    for n = 2:N-1
       Ex(n,m+1) = Ex(n,m) + dt./eps .* ( -(Hy(n,m) - Hy(n-1,m))./dz ); 
    end
    
    
end

N_mid = round(N/2);

t_plot = 2*dt:dt:(Mt+1)*dt;
t_plot = t_plot - 2.*dt;

figure(702);
subplot(2,1,1)
imagesc(t_plot, 1:N, Ex);
title('1-D FDTD Sim: E_x(z,t)', 'FontSize', 16)
xlabel('time (s)', 'FontSize', 16)
ylabel('Position on z-axis (norm units)', 'FontSize', 16)
set(gca,'Fontsize',16);
colorbar;
grid on;

subplot(2,1,2)
imagesc(t_plot, 1:N, Hy);
title('1-D FDTD Sim: H_y(z,t)', 'FontSize', 16)
xlabel('time (s)', 'FontSize', 16)
ylabel('Position on z-axis (norm units)', 'FontSize', 16)
set(gca,'Fontsize',16);
colorbar;
grid on;

figure;
plot(t_plot, Ex(N_mid,2:end))

%% Determine Frequency Spectrum - Ex

x = Ex(N_mid,2:end);
P = nextpow2(length(x));

Nq = 2.^(P+1);
Fs = 1./dt;

%% Fourier Transform:
X = fftshift(fft(x,Nq));
%% Frequency specifications:
dF = Fs/Nq;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz
%% Plot the spectrum:
figure(703);
subplot(2,2,2)
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|E_x(f)|','Fontsize',14);
title('|E_x(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 max(f)])

%figure;
subplot(2,2,1)
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|E_x(f)|','Fontsize',14);
title('|E_x(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 0.05])

figure(675);
subplot(2,1,1)
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|E_x(f)|','Fontsize',14);
title('|E_x(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 0.05])


%% Determine Frequency Spectrum - Hy

x = Hy(N_mid,2:end);
P = nextpow2(length(x));

Nq = 2.^(P+1);
Fs = 1./dt;

%% Fourier Transform:
X = fftshift(fft(x,Nq));
%% Frequency specifications:
dF = Fs/Nq;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz
%% Plot the spectrum:
figure(703);
subplot(2,2,4);
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|H_y(f)|','Fontsize',14);
title('|H_y(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 max(f)])

% figure;
subplot(2,2,3);
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|H_y(f)|','Fontsize',14);
title('|H_y(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 0.05])

figure(675);
subplot(2,1,2);
stem(f,abs(X)/Nq);
xlabel('Frequency (in hertz)','Fontsize',14);
ylabel('|H_y(f)|','Fontsize',14);
title('|H_y(f)| - Magnitude Response','Fontsize',14);
grid on;
set(gca,'Fontsize',14);
xlim([0 0.05])


%% Extra Credit 
% Use spatial FFT to extra spatial depenence of modes

Ex_spat = Ex(:, 481);

figure;
stem((1:N)*dz, Ex_spat);

% Ns = 64;
% F = fft(Ex_spat, Ns);
% 
% F2 = abs(fftshift(F));
% omega_plot = 2.*pi./32.*(0:Ns-1);
% 
% figure;
% stem(F2(1:end))
% 
% figure;
% stem(omega_plot(1:32), F2(33:end))

x = Ex(:, 481);
fs = 1/dz;

% m = length(x);          % Window length
% n = pow2(nextpow2(m));
% n = n+1;% Transform length
% y = fft(x,n);           % DFT
% f = (0:n-1)*(fs/n);     % Frequency range
% power = y.*conj(y)/n;   % Power of the DFT
% 
% y0 = fftshift(y);          % Rearrange y values
% f0 = (-n/2:n/2-1)*(fs/n);  % 0-centered frequency range
% power0 = y0.*conj(y0)/n;   % 0-centered power
% 
% figure;
% stem(f0,power0)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title('{\bf 0-Centered Periodogram}')
% %xlim([0 0.05])
% grid on;




x = Ex(:, 801);
Fs = 1/dz;


Ex_source = Ex(:, 801);

%load d_filter
%Ex_source = filter(d,Ex_source)
%Hz_source = Hz_source(571:end);
%t_plot = t_plot(571:end);

% Hz_source = cos(2.*pi.*15.*t_plot);

nfft = 2.^(2+nextpow2(length(Ex_source)));
X = fft(Ex_source,nfft);
X = X(1:nfft/2);
% Take the mag of fft of x
mx = abs(X);
% Frequency vector
f = (0:nfft/2-1)*Fs/nfft;

% Generate the plot, title and labels.


figure(462);
subplot(2,1,1);
stem(f,mx);
title('Spatial FFT Filtering of E_x(z)','Fontsize',16);
xlabel('Frequency (rad/m)','Fontsize',16);
ylabel('Mag','Fontsize',16);
xlim([0 0.1])
grid on;
set(gca,'Fontsize',16);




x = Hy(:, 801);
Fs = 1/dz;


Ex_source = Hy(:, 801);


%Ex_source = filter(d,Ex_source)

%Hz_source = Hz_source(571:end);
%t_plot = t_plot(571:end);

% Hz_source = cos(2.*pi.*15.*t_plot);

%load d_filter
%Ex_source = filter(d,Ex_source)

nfft = 2.^(2+nextpow2(length(Ex_source)));
X = fft(Ex_source,nfft);
X = X(1:nfft/2);
% Take the mag of fft of x
mx = abs(X);
% Frequency vector
f = (0:nfft/2-1)*Fs/nfft;

% Generate the plot, title and labels.


figure(462);
subplot(2,1,2);
stem(f,mx);
title('Spatial FFT Filtering of H_y(z)','Fontsize',16);
xlabel('Frequency (rad/m)','Fontsize',16);
ylabel('Mag','Fontsize',16);
xlim([0 0.1])
grid on;
set(gca,'Fontsize',16);



% Fstop1 = 0.005;
% Fpass1 = 0.01;
% Fpass2 = 0.2;
% Fstop2 = 0.21;
% Astop1 = 65;
% Apass  = 0.5;
% Astop2 = 65;
% Fs = 1/dz;
% 
% d = designfilt('bandpassfir', ...
%   'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
%   'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
%   'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
%   'StopbandAttenuation2',Astop2, ...
%   'DesignMethod','equiripple','SampleRate',Fs);
% 
% fvtool(d)