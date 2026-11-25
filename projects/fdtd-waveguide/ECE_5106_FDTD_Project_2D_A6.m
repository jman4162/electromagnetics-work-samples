%% John Hodge
%  03/13/16
% ECE 5106 FDTD Project: 1-D

clear all; clc; close all;

%% Initialize Constants

boundary = 'PEC';
eps_0 = 1;
mu_0 = 1;
c = 1/sqrt(eps_0*mu_0);
L = 30;
N = 31;
dz = 1;
h = dz;
Nx = 51; % Number of Steps in X
Ny = 51; % Number of Steps in Y
Lx = 0.1; % Length of Waveguide in X
Ly = 0.05; % Length of Waveguide in Y


Nx_mid = round(Nx/2);
Ny_mid = round(Ny/2);

dx = Lx ./ (Nx-1);
dy = Ly ./ (Ny-1);

%dt = 0.4;
n = 6;

%% Calculate Stability for FDTD
t_stab = 1 ./ (c .* sqrt(1/(dx).^2 + 1/(dy).^2));

dt = 0.4*t_stab
% dt = 3.5777e-04;

%dt = 0.3;

Mt = 10000;
t_steps = Mt;

t_plot = (0:dt:dt*(Mt-1));


%% Create Gaussian Pulse to Excite Waveguide

% clear all; clc; close all;

fs=100; %sampling frequency
sigma=1.0;
t=-0.5:1/fs:0.5; %time base
t = linspace(-pi, pi, 101);

variance=sigma^2;
x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
x_gauss_pulse = 10.*x;
subplot(2,1,1)
plot(t,x,'b');
title(['Gaussian Pulse \sigma=', num2str(sigma),'s'], 'FontSize', 16);
xlabel('Time(s)', 'FontSize', 16);
ylabel('Amplitude', 'FontSize', 16);
grid on;
set(gca,'Fontsize',16);

L=length(x);
NFFT = 1024;
X = fftshift(fft(x,NFFT));
Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
f = fs*(-NFFT/2:NFFT/2-1)/NFFT; %Frequency Vector

gauss_pulse = x;

subplot(2,1,2)
plot(f,abs(X)/(L),'r');
title('Magnitude of FFT', 'FontSize', 16);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Magnitude |X(f)|', 'FontSize', 16);
xlim([-10 10])
grid on;
set(gca,'Fontsize',16);

%% Initialize Fields

Ex = zeros(Nx,Ny,t_steps); Ey = zeros(Nx,Ny,t_steps); Ez = zeros(Nx,Ny,t_steps); Hx = zeros(Nx,Ny,t_steps); Hy = zeros(Nx,Ny,t_steps); Hz = zeros(Nx,Ny,t_steps);

% Hz(Nx_mid,Ny_mid,1) = 5;
%Hz(20,20,1:250) = 1.*ones(1,250);

%% Create Initial Mode Distribution

m = 1;
n = 0;

%Nx = 200;
%Ny = 100;

a = 0.1; % 10cm
b = 0.05; % 5cm

x_plot = linspace(0, a, Nx);
y_plot = linspace(0, b, Ny);

kx = (m*pi)/a;
ky = (n*pi)/b;
beta = 0;

w = 2*pi*4.5;
t = 0;

H0 = 1;

Hz(:,:,2) = zeros(Nx,Ny);

for xx = 2:length(x_plot)-1
    for yy = 2:length(y_plot)-1
        
        x = x_plot(xx);
        y = y_plot(yy);
        
        Hz(xx,yy,2) = H0.*cos(kx.*x).*cos(ky.*y);%.*exp(1i.*w.*t);
        
    end
end



figure;
imagesc(x_plot, y_plot, real(Hz(:,:,2))')
colorbar;
title('Initial Field Distribution', 'FontSize', 16);
xlabel('x-position (m)', 'FontSize', 16);
ylabel('y-position (m)', 'FontSize', 16);

%%


%% Step 1: Assign Initial Field Distributions for both Ex and Hy

% Define Source
% Source Location
source_x = Nx_mid;
source_y = Ny_mid;

% Main FDTD Loop
for m = 2:1:Mt
    m
    
%     if ( m < 50)
%        Hz(round(Nx_mid/2), round(Ny_mid/2)+m, 2) = 25; 
%     end
%     
%     if ( m < 50)
%        Hz(round(3*Nx_mid/2), round(3*Ny_mid/2)-m, 2) = -25; 
%     end
    
    t = m*dt; % Current time
    
    %% Step 2: Evolve EM Fields Using Finite Difference Equation
    
    % Update H-field from E-field
    for ii = 1:Nx-1
        for jj = 1:Ny-1
            
            % Hx
            Hx(ii,jj,m+1) = Hx(ii,jj,m) - dt./mu_0 .* ( (Ez(ii,jj+1,m) - Ez(ii,jj,m) )./dy + 1i.*beta.* Ey(ii,jj,m) );
            
            % Hy
            Hy(ii,jj,m+1) = Hy(ii,jj,m) + dt./mu_0 .* ( (Ez(ii+1,jj,m) - Ez(ii,jj,m) )./dx + 1i.*beta.* Ex(ii,jj,m) );
            
            % Hz
            Hz(ii,jj,m+1) = Hz(ii,jj,m) + dt./mu_0 .* ( ( Ex(ii,jj+1,m) - Ex(ii,jj,m) )./dy - ( Ey(ii+1,jj,m) - Ey(ii,jj,m) )./dx );
            
        end
    end
    
    %% Update H-field Boundaries
    
    if (strcmp(boundary, 'PEC'))
        % Upper X
        Hx(2:end-1,Ny,m) = 0; Hy(2:end-1,Ny,m) = 0; Hz(2:end-1,Ny,m) = 0;
        Ex(2:end-1,Ny,m) = 0; Ey(2:end-1,Ny,m) = 0; Ez(2:end-1,Ny,m) = 0;
        
        % Lower X
        Hy(2:end-1,1,m) = 0;
        Ex(2:end-1,1,m) = 0; Ez(2:end-1,1,m) = 0;
        
        % Left Y
        Hx(1,2:end-1,m) = 0;
        Ey(1,2:end-1,m) = 0; Ez(1,2:end-1,m) = 0;
        
        % Right Y
        Hx(Nx,2:end-1) = 0; Hy(Nx,2:end-1) = 0; Hz(Nx,2:end-1) = 0;
        Ex(Nx,2:end-1,m) = 0; Ey(Nx,2:end-1,m) = 0;  Ez(Nx,2:end-1,m) = 0;
        
        % (1,1)
        Hx(1,1,m) = 0; Hy(1,1,m) = 0;
        Ex(1,1,m) = 0; Ey(1,1,m) = 0; Ez(1,1,m) = 0;
        
        % (Nx,1)
        Hx(Nx,1,m) = 0; Hy(Nx,1,m) = 0; Hz(Nx,1,m) = 0;
        Ex(Nx,1,m) = 0; Ey(Nx,1,m) = 0; Ez(Nx,1,m) = 0;
        
        % (1,Ny)
        Hx(1,Ny,m) = 0; Hy(1,Ny,m) = 0; Hz(1,Ny,m) = 0;
        Ex(1,Ny,m) = 0; Ey(1,Ny,m) = 0; Ez(1,Ny,m) = 0;
        
        % (Nx,Ny)
        Hx(Nx,Ny,m) = 0; Hy(Nx,Ny,m) = 0; Hz(Nx,Ny,m) = 0;
        Ex(Nx,Ny,m) = 0; Ey(Nx,Ny,m) = 0; Ez(Nx,Ny,m) = 0;
        
    end
    
    %% Update E-field from H-field
    for ii = 2:Nx-1
        for jj = 2:Ny-1
            
            % Ex
            Ex(ii,jj,m+1) = Ex(ii,jj,m) + dt./eps_0 .* ( ( Hz(ii,jj,m+1) - Hz(ii,jj-1,m+1) )./dy + 1i.*beta.*Hy(ii,jj,m+1) );
            
            % Ey
            Ey(ii,jj,m+1) = Ey(ii,jj,m) - dt./eps_0 .* ( ( Hz(ii,jj,m+1) - Hz(ii-1,jj,m+1) )./dx + 1i.*beta.*Hx(ii,jj,m+1) );
            
            % Ez
            Ez(ii,jj,m+1) = Ez(ii,jj,m) + dt./eps_0 .* ( ( Hy(ii,jj,m+1) - Hy(ii-1,jj,m+1) )./dx - ( Hx(ii,jj,m+1) - Hx(ii,jj-1,m+1) )./dy );
            
        end
    end
    
    %% Update Ex Boundary
    
    if (strcmp(boundary, 'PEC'))
        % Upper X
        Hx(2:end-1,Ny,m) = 0; Hy(2:end-1,Ny,m) = 0; Hz(2:end-1,Ny,m) = 0;
        Ex(2:end-1,Ny,m) = 0; Ey(2:end-1,Ny,m) = 0; Ez(2:end-1,Ny,m) = 0;
        
        % Lower X
        Hy(2:end-1,1,m) = 0;
        Ex(2:end-1,1,m) = 0; Ez(2:end-1,1,m) = 0;
        
        % Left Y
        Hx(1,2:end-1,m) = 0;
        Ey(1,2:end-1,m) = 0; Ez(1,2:end-1,m) = 0;
        
        % Right Y
        Hx(Nx,2:end-1) = 0; Hy(Nx,2:end-1) = 0; Hz(Nx,2:end-1) = 0;
        Ex(Nx,2:end-1,m) = 0; Ey(Nx,2:end-1,m) = 0;  Ez(Nx,2:end-1,m) = 0;
        
        % (1,1)
        Hx(1,1,m) = 0; Hy(1,1,m) = 0;
        Ex(1,1,m) = 0; Ey(1,1,m) = 0; Ez(1,1,m) = 0;
        
        % (Nx,1)
        Hx(Nx,1,m) = 0; Hy(Nx,1,m) = 0; Hz(Nx,1,m) = 0;
        Ex(Nx,1,m) = 0; Ey(Nx,1,m) = 0; Ez(Nx,1,m) = 0;
        
        % (1,Ny)
        Hx(1,Ny,m) = 0; Hy(1,Ny,m) = 0; Hz(1,Ny,m) = 0;
        Ex(1,Ny,m) = 0; Ey(1,Ny,m) = 0; Ez(1,Ny,m) = 0;
        
        % (Nx,Ny)
        Hx(Nx,Ny,m) = 0; Hy(Nx,Ny,m) = 0; Hz(Nx,Ny,m) = 0;
        Ex(Nx,Ny,m) = 0; Ey(Nx,Ny,m) = 0; Ez(Nx,Ny,m) = 0;
        
    end

end

%%

% Calc Min and Max Fields
Hz_max = max(max(max(Hz)));
Hz_min = min(min(min(Hz)));

figure;
stem(squeeze(Hz(20,20,:)))

figure;
stem(squeeze(Hz(Nx_mid,Ny_mid,:)))

pts_plot = 200;

figure;
subplot(1,2,1)
imagesc(squeeze(Hz(:,Ny_mid,1:pts_plot)))
colorbar;

subplot(1,2,2)
imagesc(squeeze(Hz(Nx_mid,:,1:pts_plot)))
colorbar;



figure;
subplot(2,3,1)
imagesc(squeeze(Hz(:,Ny_mid,1:pts_plot)))
colorbar;

subplot(2,3,4)
imagesc(squeeze(Hz(Nx_mid,:,1:pts_plot)))
colorbar;

subplot(2,3,2)
imagesc(squeeze(Ex(:,Ny_mid,1:pts_plot)))
colorbar;

subplot(2,3,5)
imagesc(squeeze(Ex(Nx_mid,:,1:pts_plot)))
colorbar;

subplot(2,3,3)
imagesc(squeeze(Ey(:,Ny_mid,1:pts_plot)))
colorbar;

subplot(2,3,6)
imagesc(squeeze(Ey(Nx_mid,:,1:pts_plot)))
colorbar;

figure;
stem(squeeze(Hz(source_x,source_y,:)))
title('Hz at Source Over Time')

Hz_source_FFT = abs(fftshift(fft(squeeze(Hz(source_x,source_y,:)))));
figure;
stem(Hz_source_FFT)
%%

%% Calc and Plot FFT - Hz

Fs = 1/dt;
t_plot = (0:dt:dt*(Mt-1));

Hz_source = squeeze(Hz(source_x,source_y,:));
%Hz_source = Hz_source(571:end);
%t_plot = t_plot(571:end);

% Hz_source = cos(2.*pi.*15.*t_plot);

nfft = 2.^(2+nextpow2(length(Hz_source)));
X = fft(Hz_source,nfft);
X = X(1:nfft/2);
% Take the mag of fft of x
mx = abs(X);
% Frequency vector
f = (0:nfft/2-1)*Fs/nfft;

% Generate the plot, title and labels.
figure;
plot(t_plot,Hz_source(2:end));
title('Sine Wave Signal', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('Amplitude', 'FontSize', 16);

figure;
stem(f,mx);
title('Power Spectrum of a Sine Wave', 'FontSize', 16);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Power', 'FontSize', 16);

figure(832);
subplot(2,1,2)
stem(f,mx);
title('Hz - Power Spectrum of a Sine Wave', 'FontSize', 16);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Power', 'FontSize', 16);
grid on;
xlim([0 35]);

%% Calc and Plot FFT - Ex

Fs = 1/dt;
t_plot = (0:dt:dt*(Mt-1));

Ex_source = squeeze(Ex(source_x,source_y,:));
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
figure;
plot(t_plot,Ex_source(2:end));
title('Sine Wave Signal', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylabel('Amplitude', 'FontSize', 16);

figure;
stem(f,mx);
title('Power Spectrum of a Sine Wave');
xlabel('Frequency (Hz)');
ylabel('Power');

figure(832);
subplot(2,1,1)
stem(f,mx);
title('Ex - Power Spectrum of a Sine Wave', 'FontSize', 16);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('Power', 'FontSize', 16);
grid on;
xlim([0 35]);


%% Animate Fields

animate_fields = 0;
if (animate_fields == 1)
    for tt = 1:t_steps
        
        figure(800);
        subplot(1,3,1)
        imagesc(squeeze(Hz(:,:,tt)), [-1 1])
        title('Hz', 'FontSize', 16)
        colorbar;
        
        subplot(1,3,2)
        imagesc(squeeze(Ex(:,:,tt)), [-1 1])
        title('Ex', 'FontSize', 16)
        colorbar;
        
        subplot(1,3,3)
        imagesc(squeeze(Ey(:,:,tt)), [-1 1])
        title('Ey', 'FontSize', 16)
        colorbar;
        
    end
end

%%

% figure(180);
% colorbar;
% title('H_y Field Distribution', 'FontSize', 16);
% xlabel('x-position (m)', 'FontSize', 16);
% ylabel('y-position (m)', 'FontSize', 16);

animate_fields = 1;
if (animate_fields == 1)
    for tt = 1:t_steps
        
        figure(180);
        imagesc(x_plot, y_plot, squeeze(Hz(:,:,tt))', [-1 1])
        colorbar;
%         title('H_y Field Distribution', 'FontSize', 16);
%         xlabel('x-position (m)', 'FontSize', 16);
%         ylabel('y-position (m)', 'FontSize', 16);
%         set(gca,'Fontsize',16);
        
    end
end

% figure;
% subplot(2,1,1)
% imagesc(t_plot, 1:31, Ex)
% colorbar;
% title('Ex Field: Position vs. Time')
% ylabel('z-position')
% xlabel('time')
%
% subplot(2,1,2)
% imagesc(t_plot, 1:31, Hy)
% colorbar;
% ylabel('z-position')
% xlabel('time')
% title('Hy Field: Position vs. Time')
%
%
% figure;
% subplot(2,1,1)
% imagesc(t_plot, 1:31, Ex)
% colorbar;
% title('Ex Field: Position vs. Time')
% ylabel('z-position')
% xlabel('time')
% xlim([0 1000])
%
% subplot(2,1,2)
% imagesc(t_plot, 1:31, Hy)
% colorbar;
% ylabel('z-position')
% xlabel('time')
% title('Hy Field: Position vs. Time')
% xlim([0 1000])

% figure;
% subplot(2,1,1)
% imagesc(abs(Ex))
% colorbar;
% title('Ex Field: Position vs. Time')
% ylabel('z-position')
% xlabel('time')
%
% subplot(2,1,2)
% imagesc(abs(Hy))
% colorbar;
% ylabel('z-position')
% xlabel('time')
% title('Hy Field: Position vs. Time')
%
% figure;
% imagesc(abs(Ex).^2 + abs(Hy).^2)
% colorbar;
% title('Ex Field: Position vs. Time')
% ylabel('z-position')
% xlabel('time')


%% Animate Fields Over Time

animate_fields = 0;

if (animate_fields == 1)
    for tt = 1:t_steps,
        figure(101);
        subplot(2,1,1)
        stem(Ex(:,tt));
        axis([1 N -1 1])
        xlabel('z-position')
        ylabel('Ex Field Amplitude')
        grid on;
        drawnow
        
        subplot(2,1,2)
        stem(Hy(:,tt));
        axis([1 N -1 1])
        xlabel('z-position')
        ylabel('Hy Field Amplitude')
        grid on;
        drawnow
    end
end


% %% Step 3: Perform Fourier Transform of the EM Fields and obtain the EM field as function of frequency
%
% % Ex FFT
% Ex_FFT = fft(Ex(16,:));
% Ex_FFT = Ex_FFT./max(Ex_FFT);
%
% x = Ex(16,:);
% fs = 1./dt;
% N = length(x);
% ws = 2*pi/N;
% wnorm = -pi:ws:pi;
% wnorm = wnorm(1:length(x));
% w = wnorm*fs;
% figure;
% plot(w,abs(fftshift(Ex_FFT)))
% % axis([-30,30,0,160])
%
% % Hy FFT
% Hy_FFT = fft(Hy(16,:));
% Hy_FFT = Hy_FFT./max(Hy_FFT);
%
% figure;
% subplot(2,1,1)
% stem(abs(fftshift(Ex_FFT)))
% title('Ex: FT')
%
% subplot(2,1,2)
% stem(abs(fftshift(Hy_FFT)))
% title('Hy: FT')
%
%
% %% Processing JH 04/14/16
%
% Ex_mid = Ex(6,:);
%
% Ex_mid_FFT = fft(Ex(16,:));
% Ex_mid_FFT_norm = Ex_FFT./max(Ex_FFT);
%
% Ex_mid_FFT = fftshift(Ex_mid_FFT);
% Ex_mid_FFT_norm = fftshift(Ex_mid_FFT_norm);
%
% figure;
% plot(abs(Ex_mid_FFT_norm))
% title('Ex_mid_FFT_norm')
%
% %%
% % find out signal's Fourier transform and plot its frequency spectrum
% omega_dn = -10;
% omega_up = 10;
% N_omega = t_steps;
%
% omega_array = linspace(omega_dn,omega_up,N_omega);
% E_omega_array = zeros(1,N_omega);
%
% % Massage Variable Names
% t_array = t_plot;
% E_array = Ex_mid;
% delta_t = dt;
%
% delta_omega = omega_array(2) - omega_array(1);
% for ii=1:N_omega
%     omega = omega_array(ii);
%     exp_function_temp = exp(-1j*omega*t_array);
%     temp_vec = E_array .* exp_function_temp;
%     E_omega_array(ii) = sum(temp_vec)*delta_t;
% end
%
%
% figure;
% set(gca,'Fontsize',20);
% semilogy(omega_array/2/pi,abs(E_omega_array).^2,'.'); grid on;
% xlabel('Frequency (Hz)');
% ylabel('|E_x(\nu)|^2');
%
% figure;
% set(gca,'Fontsize',20);
% stem(omega_array/2/pi,abs(E_omega_array).^2,'.'); grid on;
% xlabel('Frequency (Hz)');
% ylabel('|E_x(\nu)|^2');