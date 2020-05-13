%% This script demonstrates chemical-shift EPI artifacts.
%
% DBE@STANFORD.EDU (May 2020) for Rad229
close all; clear all;
%% Default parameters
N=128;                    % Matrix is NxM (then padded later to accomodate motion)
gamma_bar=42.57e6;        % 1H gyromagnetic ratio [Hz/T]
B0 = 3;                   % Main (B0) field [T]
dB0_ppm = 2e-6;           % B0 shift [ppm]
ESP = 500e-6;                 % Echo spacing [s]

%% Chemical shift parameters
fat_ppm = -3.3e-6;          % Chemical shift of fat [ppm];
dB0_fat = fat_ppm*B0;       % B0 shift [T]
f_fat = gamma_bar*dB0_fat;  % Frequency offset [Hz]

%% Create an object
Obj0 = phantom('Modified Shepp-Logan',N);

%% Create a time series of phantoms with FAT phase accumulation and noise
Obj0_f = Obj0.*(Obj0==1);   % Use the brighest signals as FAT
Obj0_f = repmat(Obj0_f,[1 1 N]);
for n=1:N
  phs_f = exp(-1i*2*pi*f_fat*ESP*n);  % Phase from fat per echo
  Obj0_f(:,:,n)=Obj0_f(:,:,n).*phs_f; % Object accumulates phase during EPI duration
end

%% Create a noisy FAT object and k-space
Obj0_f = Obj0_f + (randn(size(Obj0_f))+1i*randn(size(Obj0_f)))/50; % Add some complex noise
kspc_f = zeros(size(Obj0_f));
for n=1:N
  kspc_f(:,:,n)=fftshift(fft2(Obj0_f(:,:,n))); % Define the k-space for each time
end

%% Simulate the TOP-DOWN EPI trajectory...grab one k-line from each time, sign reverse
k_epi_f = zeros(N);
for n=1:N
  k_epi_f(n,:) = kspc_f(n,:,n); % Grab one echo per ky-line during fat phase accumulation
end
IM_f = ifft2(fftshift(k_epi_f));             % Image of FAT object

%% Create a WATER object and it's k-space
Obj0_w=Obj0.*(Obj0~=1);     % Use everything but the brighest signals as WATER
Obj_w = Obj0_w + (randn(size(Obj0_w))+i*randn(size(Obj0_w)))/25; % Add some complex noise

kspc_w=fftshift(fft2(Obj_w));   % Define the k-space
IM_w = ifft2(fftshift(kspc_w)); % Image of WATER object

%% The FT is a linear operator so we can add the k-space or the images
IM_CS = IM_w+IM_f;

%% Create a figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0));     title('Object');  colorbar;
subplot(2,2,2); imagesc(abs(IM_CS));    title('Image');   colorbar;
subplot(2,2,3); imagesc(abs(k_epi_f));  title('Phase');   colorbar;
subplot(2,2,4); imagesc(angle(IM_CS));  title('Fat k-space'); colorbar;