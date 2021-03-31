%% Rad229_Fourier_Encoding demonstrates some Fourier sampling functions
%
% SYNTAX - F = Rad229_Fourier_Encoding(acq, n_kx, n_ky, upsamp);
%
% INPUTS -   acq.FOVx = 210e-3;  % Field-of-view along y-direction [m]
%            acq.Nx = 21;        % Number of pixels to discretize FOVy [#]
%            acq.RBW = 32e3;     % Receiver pixel bandwidth magnitude [Hz] 
% 
%            acq.FOVy = 150e-3;  % Field-of-view along y-direction [m]
%            acq.Ny = 15;        % Number of pixels to discretize FOVy [#]
% 
%            n_kx = 0;           % Fourier kx-sampling point to show
%            n_ky = 4;           % Fourier ky-sampling point to show
%
%            upsamp = 5;        % Upsampling factor for visualizing patterns
%
% OUTPUTS -  F - Fourier encoding pattern
%
% DBE@STANFORD.EDU (March 2021) for Rad229

%% These questions can be used to further explore the code and concepts.
%
% 1) Set n_kx = 3 and n_ky = 0. How many phase cycles do you see along the frequency encode direction? 
%
% 2) Set n_kx = 0 and n_ky = 4. How many phase cycles do you see along the phase encode direction? 
%
% 3) Set n_kx = 3 and n_ky = 5. How many phase cycles do you see along the phase encode direction? 
%
% 4) Simplify the code below to skip the gradient amplitude calculation.
%
% 5) [Advanced] Use this code to estimate the Fourier coefficients for an object.

function F = Rad229_Fourier_Encoding(n_kx, n_ky, upsamp)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Define acquisition parameters
if nargin == 0 
  acq.FOVx = 210e-3;  % Field-of-view along y-direction [m]
  acq.Nx = 21;        % Number of pixels to discretize FOVy [#]
  acq.RBW = 32e3;     % Receiver pixel bandwidth magnitude [Hz] 

  acq.FOVy = 150e-3;  % Field-of-view along y-direction [m]
  acq.Ny = 15;        % Number of pixels to discretize FOVy [#]
    
  n_kx = 3;
  n_ky = 5;

  upsamp = 5;        % Upsampling factor for visualizing patterns
end

%% Calculate some spatial parameters
acq.dx = acq.FOVx / acq.Nx;     % Pixel dimension along y-direction [m]
acq.x_pos = linspace(-acq.FOVx/2, acq.FOVx/2, upsamp * acq.Nx); % Define the pixel locations [m]
acq.dy = acq.FOVy / acq.Ny;     % Pixel dimension along y-direction [m]
acq.y_pos = linspace(-acq.FOVy/2, acq.FOVy/2, upsamp * acq.Ny); % Define the pixel locations [m]

[acq.X, acq.Y] = ndgrid(acq.x_pos, acq.y_pos); % Define a grid of pixel positions [m]

%% Compute the applied gradient
Gx.G_plat_amp = ( 2 * acq.RBW ) / ( sys.gamma_bar * acq.FOVx );  % Readout gradient plateau amplitude [T/m]
Gy.G_plat_amp = ( 2 * acq.RBW ) / ( sys.gamma_bar * acq.FOVy );  % Readout gradient plateau amplitude [T/m]
  % Note - This assume the phase encode step takes as long as a readout step

%% Calculate the delta k-space steps
acq.dkx = sys.gamma_bar * Gx.G_plat_amp * ( 1 / (2 * acq.RBW) );       % kx-space component [1/m]
acq.dky = sys.gamma_bar * Gy.G_plat_amp * ( 1 / (2 * acq.RBW) );       % ky-space component [1/m]

%% Compute the Fourier sampling functions
F = exp( -1i * 2 * pi * (n_kx * acq.dkx * acq.X + n_ky * acq.dky * acq.Y) ); % Fourier sampling functions

%% Display the sampling function
f = figure; hold on;
  s = surf(acq.X, acq.Y, angle(F)); 
  view(0,90); set(s,'EdgeColor','None'); colorbar; caxis([-pi pi]);
  title('angle(F)'); axis image xy; set(gca, 'FontSize', 24);