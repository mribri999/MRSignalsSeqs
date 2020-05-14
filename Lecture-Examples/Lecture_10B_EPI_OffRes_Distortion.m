%% This script demonstrates off-resonance distortion in EPI.
%
% DBE@STANFORD.EDU (May 2020) for Rad229
close all; clear all;
%% Default parameters
N=128;                    % Matrix is NxM (then padded later to accomodate motion)
gamma_bar=42.57e6;        % 1H gyromagnetic ratio [Hz/T]
ESP = 500e-6;                 % Echo spacing [s]

%% Create a noisy object and k-space
Obj0 = phantom('Modified Shepp-Logan',N);

%% Create a time series of phantoms with off-resonance phase accumulation and noise
Obj0_dist = repmat(Obj0,[1 1 N]);
freq=linspace(0,5,N).^2;  % Arbitrary off-resonance frequency [Hz]
freq_scl = freq'*freq;    % Define an off-resonance frequency across the FOV

for n=1:N
  phs_off = exp(-1i*2*pi*freq_scl*ESP*n);     % Phase from off-resonance per echo
  Obj0_dist(:,:,n)=Obj0_dist(:,:,n).*phs_off; % Object accumulates phase during EPI duration
end

% %%
% % disp=2*pi*freq_scl*ESP;
% disp=freq_scl*ESP*N;
% figure; imagesc(disp); colorbar;

%% Create a noisy object and k-space
Obj0_dist = Obj0_dist + (randn(size(Obj0_dist))+1i*randn(size(Obj0_dist)))/50; % Add some complex noise
kspc_dist = zeros(size(Obj0_dist));
for n=1:N
  kspc_dist(:,:,n)=fftshift(fft2(Obj0_dist(:,:,n))); % Define the k-space for each time
end

%% Simulate the TOP-DOWN EPI trajectory...grab one k-line from each time
k_epi_f = zeros(N);
for n=1:N
  kspc_epi_dist(n,:) = kspc_dist(n,:,n); % Grab one echo per ky-line during phase accumulation
end
IM_epi_dist = ifft2(fftshift(kspc_epi_dist)); % Image of distorted object

%% Create a figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0)); title('Object'); colorbar;
subplot(2,2,2); imagesc(abs(IM_epi_dist)); title('Image'); colorbar;
subplot(2,2,3); imagesc(freq_scl); title('Off-Resonance [Hz]'); colorbar;
subplot(2,2,4); imagesc(angle(IM_epi_dist)); title('Off-Resonance [rad]'); colorbar;

%% Simulate the BOTTOM-UP EPI trajectory...grab one k-line from each time
k_epi_f = zeros(N);
for n=1:N
  kspc_epi_dist(N-(n-1),:) = kspc_dist(N-(n-1),:,n); % Fill the bottom ky-line (N-(n-1)) with the first echo (n)
end
IM_epi_dist2 = ifft2(fftshift(kspc_epi_dist)); % Image of distorted object

figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0)); title('Object'); colorbar;
subplot(2,2,2); imagesc(abs(IM_epi_dist)); title('Image (Top-Down)'); colorbar;
subplot(2,2,3); imagesc(freq_scl); title('Off-Resonance [Hz]'); colorbar;
subplot(2,2,4); imagesc(abs(IM_epi_dist2)); title('Image (Bottom-Up)'); colorbar;
