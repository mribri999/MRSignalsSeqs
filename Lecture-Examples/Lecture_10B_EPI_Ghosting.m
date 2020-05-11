%% This script demonstrates ghosting artifacts in EPI.
%
% DBE@STANFORD.EDU (May 2020) for Rad229
%

%% EPI Nyquist ghosting (CONSTANT phase shift)
close all; clear all;

% Default parameters
N=128;                    % Matrix is NxM (then padded later to accomodate motion)
gamma_bar=42.57e6;        % 1H gyromagnetic ratio [Hz/T]
B0 = 3;                   % Main (B0) field [T]
dB0_ppm = 3.5e-6;           % B0 shift [ppm]
ESP = 1000e-6;                 % Echo spacing [s]
 
%% Create a noisy object
Obj0 = phantom('Modified Shepp-Logan',N/2);
Obj0 = padarray(Obj0,[N/4 N/4],0,'both'); % Just helps see aliasing...
Obj = Obj0 + (randn(size(Obj0))+i*randn(size(Obj0)))/25; % Add some complex noise

%% EPI Nyquist ghosting (LINEAR phase shift)
% Calculate a constant phase that arises from a B0-field (B0)
dB0 = dB0_ppm*B0;         % B0 shift [T]
f0 = gamma_bar*dB0;       % Frequency offset [Hz]
dphs = 2*pi*f0*ESP;        % Phase offset [rad]
phs = dphs*ones(size(Obj)); % Constant phase offset [rad]

phs_eff=zeros(size(phs));
phs_eff(2:2:end,:) = +phs(2:2:end,:);
phs_eff(1:2:end,:) = -phs(1:2:end,:);

% Even echoes "see" the +phs; odd echoes "see" the -phs
Obj_even = Obj0.*exp(-1i*phs);
Obj_odd = Obj0.*exp(+1i*phs);

%% EPI readout mixes even and odd echoes...
k_spc_even=fftshift(fft2(Obj_even)); % Define the k-space
k_spc_odd=fftshift(fft2(Obj_odd));   % Define the k-space

k_spc = zeros(size(k_spc_even));
k_spc(2:2:end,:) = k_spc_even(2:2:end,:);
k_spc(1:2:end,:) = k_spc_odd(1:2:end,:);

IM = ifft2(fftshift(k_spc)); % Image of object

%% Create a figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0));  title('Object');  colorbar; axis equal tight;
subplot(2,2,2); imagesc(abs(IM));    title('Image');   colorbar; axis equal tight;
subplot(2,2,3); imagesc(phs_eff);    title('Object Phase');   colorbar; axis equal tight;
subplot(2,2,4); imagesc(abs(k_spc)); title('k-space'); colorbar; axis equal tight; caxis([0 150]);


%% EPI Nyquist ghosting (LINEAR phase shift)
% Calculate a constant phase that arises from a B0-field (B0)
dB0 = dB0_ppm*B0;         % B0 shift [T]
f0 = gamma_bar*dB0;       % Frequency offset [Hz]
dphs = 2*pi*f0*ESP;        % Phase offset [rad]
phs = repmat(dphs*linspace(0,1,N)',[1 N]); % Constant phase offset [rad]

phs_eff=zeros(size(phs));
phs_eff(2:2:end,:) = +phs(2:2:end,:);
phs_eff(1:2:end,:) = -phs(1:2:end,:);

% Even echoes "see" the +phs; odd echoes "see" the -phs
Obj_even = Obj0.*exp(-1i*phs);
Obj_odd = Obj0.*exp(+1i*phs);

%% EPI readout mixes even and odd echoes...
k_spc_even=fftshift(fft2(Obj_even)); % Define the k-space
k_spc_odd=fftshift(fft2(Obj_odd));   % Define the k-space

k_spc = zeros(size(k_spc_even));
k_spc(2:2:end,:) = k_spc_even(2:2:end,:);
k_spc(1:2:end,:) = k_spc_odd(1:2:end,:);

IM = ifft2(fftshift(k_spc)); % Image of object

%% Create a figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0));  title('Object');  colorbar; axis equal tight;
subplot(2,2,2); imagesc(abs(IM));    title('Image');   colorbar; axis equal tight;
subplot(2,2,3); imagesc(phs_eff);    title('Object Phase');   colorbar; axis equal tight;
subplot(2,2,4); imagesc(abs(k_spc)); title('k-space'); colorbar; axis equal tight; caxis([0 150]);