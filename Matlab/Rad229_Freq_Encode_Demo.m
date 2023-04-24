%% Rad229_Freq_Encode_Demo demonstrates designing the frequency encoding gradient for imaging.
%
% In general, the user can usually pick RBW, FOVx, and Nx on the scanner. This permits the 
% calculation of the readout gradient amplitude (Gx) and duration (T_acq).
%
%    Gx = 2|RBW| / ( gamma_bar * FOVx )
%    T_acq = Nx / 2|RBW|
%
% Note: See code for variable definitions and units
%
% SYNTAX -  [acq, sys, Gx, Gy, Gz, RF] = Rad229_Freq_Encode_Demo(acq)
%
% INPUTS -  acq.RBW = 32e3;     % Receiver pixel bandwidth magnitude [Hz] 
%           acq.FOVx = 256e-3;  % Field-of-view along x-direction [m]
%           acq.Nx = 128;       % Number of pixels to discretize FOVx [#]
%
% OUTPUTS - acq - A structure of the defined acquisition parameters.
%           sys - A structure of the defined MRI system parameters.
%           Gx - A structure with the Gx gradient parameters.
%           Gy - A structure with the Gy gradient parameters.
%           Gz - A structure with the Gz gradient parameters.
%           RF - A structure with the RF gradient parameters.   
%  
%    RBW = Receiver bandwidth magnitude [Hz] 
%    Gx = Readout gradient amplitude
%    FOVx = field-of-view along x-direction [m]
%
% DBE@STANFORD.EDU (March 2021) for Rad229

%% Frequency Encode Questions:
%
%     1) How does the readout gradient change if the RBW is doubled? Halved?
%       
%     2) Increase the RBW to ±125kHz. What happens? Why? What hardware change would eliminate this problem?
%  
%     3) Decrease the FOV to 128mm. What happens?
%
%     4) [Advanced] Demonstrate k-space traversal during the readout pre-phasing and readout lobes.
%
%     5) [Advanced] Incorporate a FOV shift in the frequency encoding direction. Hint: See Bernstein p. 255

function [acq, sys, Gx, Gy, Gz, RF] = Rad229_Freq_Encode_Demo(acq)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Define Frequency Encoding Variables - Most scanners are configured for the user to select RBW, 
%  FOVx, and Nx on the scanner, which then defines the spatial resolution.
if nargin == 0
  acq.RBW = 128e3;     % Receiver pixel bandwidth magnitude [Hz] 
% **NOTE** acq.RBW is ±RBW on the scanner. A 2x-factor is used in the equations.
  acq.FOVx = 256e-3;  % Field-of-view along x-direction [m]
  acq.Nx = 128;       % Number of pixels to discretize FOVx [#]
  acq.dx = acq.FOVx / acq.Nx;     % Pixel dimension along x-direction [m]
  acq.x_pos = linspace(-acq.FOVx/2, acq.FOVx/2, acq.Nx); % Define the pixel locations [m]
end

%% Design the frequency encode gradients
Gx.G_plat_amp = ( 2 * acq.RBW ) / ( sys.gamma_bar * acq.FOVx );  % Readout gradient plateau amplitude [T/m]
  if Gx.G_plat_amp > sys.G_max, warning('Gx.G_plat_amp exceeds sys.G_max!'); end
Gx.t_acq = acq.Nx / ( 2 * acq.RBW ); % Readout gradient acquisition (plateau) duration [s]

Gx.t_ramp = sys.G_max / sys.S_max;  % Ramp time [s]
%   if ~isinteger(Gx.t_ramp / sys.dt), warning('Gx.t_ramp is NOT an integer number of sys.dt!'); end

tmp = 0 : sys.dt : Gx.t_acq;  % Readout gradient time vector [s]
Gx.dG = sys.S_max * sys.dt;  % Readout gradient incremental gradient step [T/m]
Gx.G_ramp = ( 0 : Gx.dG : Gx.G_plat_amp )'; % Readout gradient ramp [T/m]
Gx.G_plat = ones( size(tmp') ) * Gx.G_plat_amp;  % Phase encode plateau waveforms [T/m]

% Composite the readout gradient (first ramp, plateau, and last ramp)
Gx.G = [Gx.G_ramp; Gx.G_plat; flipud(Gx.G_ramp)];

%% Design the readout pre-phasing gradient
A = sum( Gx.G * sys.dt ); % The pre-phasing gradient needs to have half this area [T*s/m]
Gx.G_pre = -sys.G_max; % The pre-phasing gradient should use the maximum hardware to be fast [T/m] 
Gx.t_pre = ( A / 2 ) / abs( Gx.G_pre ) - Gx.t_ramp; % Duration of pre-phasing gradient plateau [s]

tmp = 0 : sys.dt : Gx.t_pre;  % Pre-phasing gradient gradient time vector [s]

Gx.G_pre_ramp = ( 0 : -Gx.dG : Gx.G_pre )';     % Pre-phasing gradient ramp waveform [T/m]
Gx.G_pre_plat = ones( size(tmp') ) * Gx.G_pre;  % Pre-phasing gradient plateau waveforms [T/m]

%% Composite the entire readout gradient
Gx.G = [Gx.G_pre_ramp; Gx.G_pre_plat; flipud(Gx.G_pre_ramp); Gx.G_ramp; Gx.G_plat; flipud(Gx.G_ramp)];
RF.B1 = zeros( size(Gx.G , 1) , 1 );
Gy.G  = zeros( size(Gx.G , 1) , 1 );
Gz.G  = zeros( size(Gx.G , 1) , 1 );

%% Plot the frequency encoding waveform
Rad229_PSD_fig(1e6*RF.B1, Gx.G, Gy.G, Gz.G, sys.dt);

return