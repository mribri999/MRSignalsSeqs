%% Rad229_MRI_Signal_Eqn -- INCOMPLETE FUNCTION
%
% This script demonstrates the MRI signal equation and Fourier encoding for
% an imaging experiment
%
% DBE@STANFORD.EDU (March 2021) for Rad229

%% MRI Signal Equation Questions:
%     X) Make a plot of the k-space trajectories from the gradient
%     waveforms. Distinguish between k-points that are acquired and those
%     that are not.

error('This function is incomplete.')

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Load some MRI data (default with MATLAB installation)
load MRI; D = double(D(:,:,1,10)); % Grab a single slice and convert to DOUBLE
I = imresize(D, 0.3);             % It will be helpful to use a low-res object
I = I(:,2:end-1);                  % Assymetric matrix dimensions are helpful for debugging!

%% Define Phase Encoding Variables - Most scanners are configured for the user to select the FOV and 
%  matrix size for encoding, which then defines the spatial resolution.
acq.RBW = 32e3;     % Receiver pixel bandwidth magnitude [Hz] 
% acq.RBW = 50e3;     % Receiver pixel bandwidth magnitude [Hz] 
acq.FOVx = 100e-3;  % Field-of-view along y-direction [m]
acq.Nx = size(I,2);        % Number of pixels to discretize FOVy [#]
acq.dx = acq.FOVx / acq.Nx;     % Pixel dimension along y-direction [m]
acq.x_pos = linspace(-acq.FOVx/2, acq.FOVx/2, acq.Nx); % Define the pixel locations [m]

acq.FOVy = 100e-3;  % Field-of-view along y-direction [m]
acq.Ny = size(I,1);        % Number of pixels to discretize FOVy [#]
acq.dy = acq.FOVy / acq.Ny;     % Pixel dimension along y-direction [m]
acq.y_pos = linspace(-acq.FOVy/2, acq.FOVy/2, acq.Ny); % Define the pixel locations [m]

[acq.X, acq.Y] = ndgrid(acq.x_pos, acq.y_pos); % Define a grid of pixel positions [m]

%% Design the frequency and phase encoding gradients
[acq, sys, Gx_freq, Gy_freq, Gz_freq, RF_freq] = Rad229_Freq_Encode_Demo(acq);
[acq, sys, Gx_phs , Gy_phs , Gz_phs , RF_phs ] = Rad229_Phase_Encode_Demo(acq);

%% Combine the gradient axes
Gx = [Gx_phs.G; Gx_freq.G];
Gy = [Gy_phs.G; repmat(Gy_freq.G, 1, size(Gy_phs.G, 2))];
Gz = [Gz_phs.G; Gz_freq.G];
RF = [RF_phs.B1; RF_freq.B1];

%% Plot the gradient waveforms
Rad229_PSD_fig(1e6*RF, Gx, Gy, Gz, sys.dt);

%% Calculate the k-space trajectories - Note this is the trajectory using sys.dt raster points
%  This sampling is different than the receiver systems bandwidth (RBW) sampling
kx = sys.gamma_bar * cumsum(Gx,1) * sys.dt;
ky = sys.gamma_bar * cumsum(Gy,1) * sys.dt;
kz = sys.gamma_bar * cumsum(Gz,1) * sys.dt;

%% Define the k-space points sampled using the sys.dt clock and the acq.RBW clock
ind0 = length( [Gx_phs.G; Gx_freq.G_pre_ramp; Gx_freq.G_pre_plat; Gx_freq.G_pre_ramp; Gx_freq.G_ramp] ) + 1;
ind1 = ind0 + length(Gx_freq.G_plat) - 1;

% k-points on the sys.dt clock
acq.kx = kx(ind0:ind1);
acq.ky = ky(ind0:ind1,:);
acq.kz = kz(ind0:ind1);

% k-points on the acq.RBW clock
acq.kx_rs = linspace( kx(ind0), kx(ind1), acq.Nx );
acq.ky_rs = acq.ky(:,1:acq.Nx);
acq.kz_rs = linspace( kz(ind0), kz(ind1), acq.Nx );

figure; hold on;
  plot(kx,ky,'.-');
  plot(acq.kx, acq.ky, 'ko');
  plot(acq.kx_rs, acq.ky_rs, 'g.');

%% Calculate the Fourier encoding 
upsamp = 4;
acq.x_pos_hi = linspace(-acq.FOVx/2, acq.FOVx/2, upsamp * acq.Nx); % Define the pixel locations [m]
acq.y_pos_hi = linspace(-acq.FOVy/2, acq.FOVy/2, upsamp * acq.Ny); % Define the pixel locations [m]

[acq.X_hi, acq.Y_hi] = meshgrid(acq.x_pos_hi, acq.y_pos_hi); % Define a plaid grid of pixel positions [m]

warning('Trying to hit kx=0 and ky=0 with assymetric and odd Image dimensions')

% F = zeros(acq.Nx, acq.Ny, size(acq.X_hi,1), size(acq.Y_hi,2)); % Initialize the Fourier arrays
% for nx = 1 : acq.Nx
%   for ny = 1 : acq.Ny
%     F(nx, ny, :, :) = exp( -1i * 2 * pi * (acq.kx_rs(nx).*acq.X_hi + acq.ky_rs(nx,ny).*acq.Y_hi) ); % Fourier sampling functions
%   end
% end

F = zeros(acq.Ny, acq.Nx, size(acq.Y_hi,1), size(acq.X_hi,2)); % Initialize the Fourier arrays
for nx = 1 : acq.Nx
  for ny = 1 : acq.Ny
    F(ny, nx, :, :) = exp( -1i * 2 * pi * (acq.kx_rs(nx).*acq.X_hi + acq.ky_rs(ny,nx).*acq.Y_hi) ); % Fourier sampling functions
  end
end

%% Make a figure of the Fourier encoding patterns

warning('Not updating with different 16s');

figure; hold on;
  imagescn( angle(squeeze( F(19, 18, :, :) ) ) );
%   imagescn( angle(squeeze( F(round((acq.Ny-1)/2), round((acq.Nx-1)/2), :, :) ) ) );

figure; hold on;
  imagescn( angle(squeeze( F(37, 17, :, :) ) ) );

figure; hold on;
  imagescn( angle(squeeze( F(37, 39, :, :) ) ) );

%% Demonstrate the the Fourier encoding patterns can reproduce the object

% % % dx=1;
% % % dy=1;
% % % FOVx=dx*Nx;     % Field of view along x-direction
% % % FOVy=dy*Ny;     % Field of view along y-direction
% % % dkx=1/FOVx;         % Incremental spatial frequency step along x-direction
% % % dky=1/FOVy;         % Incremental spatial frequency step along y-direction
% % 
% % % dkx=gamma_bar*Gx*dt;
% % % dky=gamma_bar*Gy*dt;
% % 
% % kx=sys.gamma_bar*sum(Gx(1:t_ind))*dt;
% % ky=sys.gamma_bar*sum(Gy(1:t_ind))*dt;
% % 
% % F=exp(-1i*2*pi*(kx*X+ky*Y)); % Fourier sampling functions
% % 
% % %  Need spatial extent --> X and Y
% % %  Need gradient history or cummulative time...probably integral form of
% % %  equation
% %     
% % % for ny=(-Ny/8):2:(Ny/8-1)   % Don't go to the edges; skip some points...
% % %   for nx=(-Nx/8):2:(Nx/8-1)  % Don't go to the edges; skip some points...
% % %     F(:,:,indx,indy)=exp(-1i*2*pi*(nx*dkx*X+ny*dky*Y)); % Fourier sampling functions
% % %     indx=indx+1;
% % %     ind=ind+1;
% % %   end
% % %   indx=1;
% % %   indy=indy+1
% % % end