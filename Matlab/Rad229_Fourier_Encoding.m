%% Rad229_Fourier_Encoding demonstrates some Fourier sampling functions
%
% SYNTAX - [acq, F] = Rad229_Fourier_Encoding(acq);
%
% INPUTS -   acq.FOVx = 210e-3;  % Field-of-view along y-direction [m]
%            acq.Nx = 21;        % Number of pixels to discretize FOVy [#]
% 
%            acq.FOVy = 150e-3;  % Field-of-view along y-direction [m]
%            acq.Ny = 15;        % Number of pixels to discretize FOVy [#]
% 
%            acq.n_kx = 0;       % Fourier kx-sampling point to show (=0 is middle of k-space)
%    
%            acq.n_ky = 4;       % Fourier ky-sampling point to show (=0 is middle of k-space)
%
%            acq.upsamp = 5;     % Upsampling factor for visualizing patterns
%
% OUTPUTS -  acq - Acquisition parameter structure
%              F - Fourier encoding pattern
%
% DBE@STANFORD.EDU (March 2021) for Rad229

%% These questions can be used to further explore the code and concepts.
%
%   1) Set acq.n_kx = 3 and acq.n_ky = 0. How many phase cycles do you see along the frequency encode direction? 
%
%   2) Set acq.n_kx = 0 and acq.n_ky = 4. How many phase cycles do you see along the phase encode direction? 
%
%   3) Set acq.n_kx = 2 and acq.n_ky = 5. How many phase cycles do you see along the phase and frequency encode direction? 
%
%   4) If you set acq.n_kx = k_x,max and acq.n_ky = k_y,max does this satisfy Nyquist sampling? Explain.
%
%   4) [Advanced] Revise the code to compute the applied phase and frequency encoding gradients, then dkx and dky.
%
%   5) [Advanced] Use this code to estimate the Fourier coefficients for an object, then use the FFT 
%      to recover an image of the object.
%
%   6) [Advanced] Use this code to demonstrate field-of-view aliasing and compressed-sensing artifacts.

function [acq, F] = Rad229_Fourier_Encoding(acq)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Define acquisition parameters
if nargin == 0 
  acq.FOVx = 210e-3;  % Field-of-view along y-direction [m]
  acq.Nx = 21;        % Number of pixels to discretize FOVy [#]
%   acq.RBW = 32e3;     % Receiver pixel bandwidth magnitude [Hz] 

  acq.FOVy = 150e-3;  % Field-of-view along y-direction [m]
  acq.Ny = 15;        % Number of pixels to discretize FOVy [#]
    
  acq.n_kx = 0;
  acq.n_ky = 11;

  acq.upsamp = 5;        % Upsampling factor for visualizing patterns
end

%% Check the index is within range
if abs(acq.n_kx) > ( (acq.Nx - 1) / 2 ), warning('acq.n_kx will exceed kx_max.'); end
if abs(acq.n_ky) > ( (acq.Ny - 1) / 2 ), warning('acq.n_ky will exceed ky_max.'); end

%% Calculate some spatial parameters
acq.x_pos = linspace(-acq.FOVx/2, acq.FOVx/2, acq.upsamp * acq.Nx); % Define the pixel locations [m]
acq.y_pos = linspace(-acq.FOVy/2, acq.FOVy/2, acq.upsamp * acq.Ny); % Define the pixel locations [m]

[acq.X, acq.Y] = ndgrid(acq.x_pos, acq.y_pos); % Define a grid of pixel positions [m]

%% Calculate the delta k-space steps
acq.dkx = 1 / acq.FOVx;                    % Spatial frequency increment [cycle/m]
acq.dky = 1 / acq.FOVy;                    % Spatial frequency increment [cycle/m]

%% Compute the Fourier sampling functions
F = exp( -1i * 2 * pi * (acq.n_kx * acq.dkx * acq.X + acq.n_ky * acq.dky * acq.Y) ); % Fourier sampling functions

%% Display the sampling function
f = figure; hold on;
  s = surf(acq.X, acq.Y, angle(F)); 
  view(0,90); set(s,'EdgeColor','None'); colorbar; caxis([-pi pi]);
  title('angle(F)'); axis image xy; set(gca, 'FontSize', 24);
  
%% Compute the applied gradient
% acq.dx = acq.FOVx / acq.Nx;     % Pixel dimension along y-direction [m]
% acq.dy = acq.FOVy / acq.Ny;     % Pixel dimension along y-direction [m]

% Gx.G_plat_amp = ( 2 * acq.RBW ) / ( sys.gamma_bar * acq.FOVx );  % Readout gradient plateau amplitude [T/m]
% Gy.G_plat_amp = ( 2 * acq.RBW ) / ( sys.gamma_bar * acq.FOVy );  % Readout gradient plateau amplitude [T/m]
  % Note - This assume the phase encode step takes as long as a readout step

% Compute dkx and dky from the applied gradients
% acq.dkx = sys.gamma_bar * Gx.G_plat_amp * ( 1 / (2 * acq.RBW) );       % kx-space component [1/m]
% acq.dky = sys.gamma_bar * Gy.G_plat_amp * ( 1 / (2 * acq.RBW) );       % ky-space component [1/m]