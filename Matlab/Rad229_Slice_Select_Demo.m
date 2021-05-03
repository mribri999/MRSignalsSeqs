
%% Rad229_Slice_Select_Demo – Demonstrate designing both the RF pulse and gradient for slice selection.
% 
%  SYNTAX - [acq, sys, Gx, Gy, Gz, RF] = Rad229_Slice_Select_Demo(acq, RF);
%
%  INPUTS -  acq.dz = 10.0e-3;         % Slice thickness [m]
%            acq.z0 = 0.0;             % Slice position [m] (relative to isocenter)
%
%            RF.alpha = 90;            % RF flip angle [degrees]
%            RF.TBW = 6;               % Time*bandwidth product [unitless] (relates to quality of RF pulse)
%            RF.dur = 2.5e-3;          % RF pulse duration [seconds]
%            RF.apod = 0.46;           % apod = 0 is non-apodized, apod = 0.5 is 'Hanning windowed', apod<0.5 is 'Hamming windowed'
%
%  OUTPUTS – acq - A structure of the defined acquisition parameters.
%            sys - A structure of the defined MRI system parameters.
%            Gx - A structure with the Gx gradient parameters.
%            Gy - A structure with the Gy gradient parameters.
%            Gz - A structure with the Gz gradient parameters.
%            RF - A structure with the RF gradient parameters.   
% 
% DBE@STANFORD.EDU (March 2021) for Rad229

%% These questions can be used to further explore the code and concepts.
% RF Questions:
%     1.1) What is the highest possible flip angle without a B1 warning?
%     
%     1.2) A 180' pulse exceeds B1,max. What adjustments can be made to enable a 180'?
%
%     1.3) Design the shortest possible TBW=8 refocusing pulse that doesn't exceed B1,max.
%
%     1.4) [Advanced] Update the code to account for the slice position.
%
%     1.5) [Advanced] ake a plot of peak B1 for the shortest duration RF pulses for flip angles 
%          from 5 to 180 degrees. All pulses should be within the default hardware limits. What do 
%          you observe about the shape of the curve and the hardware limits encountered.
%
% Gradient Questions:
%     2.1) What is the thinnest slice you can excite without a gradient warning?
%
%     2.2) What adjustments can be made to enable exciting a 0.1mm slice?
%
%     2.3) [Advanced] Demonstrate the effectiveness of the post-excitation refocusing gradient.
%
%     2.4) [Advanced] What is the thinnest slice you can excite with a 10' flip angle 
%          2ms TBW=6 RF pulse (within hardware limits)? What is the shortest duration you can 
%          make this pulse?
%
% Non-1H Imaging Questions (sys.gamma_bar=10.7084e6 [Hz/T]):
%     3.1) What happens if the default pulse is used to excite 13C? Why?
%
%     3.2) What else needs to change about the RF pulse to excite 13C?
%
%     3.3) [Advanced] Design an TBW=10 RF pulse to refocus 13C for a 10mm slice within hardware 
%          limits. If you were to build an MRI system for 13C imaging would you prefer an increase
%          in G_max or B1_max? Why?

function [acq, sys, Gx, Gy, Gz, RF] = Rad229_Slice_Select_Demo(acq, RF)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;
% sys.gamma_bar=10.7084e6;  % Gyromagnetic ratio for 11C

%% Define the RF pulse to excite a slice
if nargin == 0
  % Design parameters
  acq.dz = 10.0e-3;         % Slice thickness [m]
  acq.z0 = 0.0;             % Slice position [m] (relative to isocenter)

  RF.alpha = 90;            % RF flip angle [degrees]
  RF.TBW = 6;               % Time*bandwidth product [unitless] (relates to quality of RF pulse)
  RF.dur = 2.0e-3;          % RF pulse duration [seconds]
  RF.apod = 0.46;           % apod = 0 is non-apodized, apod = 0.5 is 'Hanning windowed', apod<0.5 is 'Hamming windowed'
end

% Calculated parameters
RF.BW = RF.TBW ./ RF.dur;            % RF bandwidth [Hz]
RF.N = ceil( RF.dur / sys.dt );      % Number of points in the pulse
RF.t = linspace( -1, 1, RF.N )';      % Normalized time vector (duration is 1)

%% Design the RF pulse (waveform)
RF.apod_fxn = ( 1 - RF.apod ) + RF.apod * cos( pi * RF.t );  % RF pulse apodization function
RF.B1 = sinc( RF.t * RF.TBW / 2 ) .* RF.apod_fxn;            % Apodized (windowed) RF SINC pulse
RF.B1 = RF.B1 ./ sum( RF.B1 );                               % Normalize to unit area
RF.B1 = ( RF.alpha * pi/180 ) * RF.B1 / ( 2 * pi * sys.gamma_bar * sys.dt); % Scale for flip angle and real amplitude
  if max( RF.B1 ) > sys.B1max, warning('RF.B1 exceeds sys.B1max'); end

%% Design the gradient to excite a slice
Gz.G_amp = RF.BW / ( sys.gamma_bar * acq.dz ); % Gradient amplitude [T/m]
Gz.G = Gz.G_amp * ones( size( RF.B1 ) );       % Gradient amplitude [T/m]
  if max(Gz.G)>sys.G_max, warning('Gz.G exceeds sys.G_max'); end

Gz.dG = sys.S_max * sys.dt;  % Slice-select gradient incremental gradient step [T/m]

Gz.G_ramp = ( 0 : Gz.dG : Gz.G_amp )'; % Readout gradient ramp [T/m]
Gz.t_ramp = sys.G_max / sys.S_max;  % Ramp time [s]
  if ~isinteger(Gz.t_ramp / sys.dt), warning('Gz.t_ramp is NOT an integer number of sys.dt!'); end

%% Design the post-excitation refocusing gradient
A = sum( Gz.G * sys.dt );  % The refocusing gradient needs to have half this area [T*s/m]
Gz.G_post = -sys.G_max;    % The refocusing gradient should use the maximum hardware to be fast [T/m] 
Gz.t_post = ( A / 2 ) / abs( Gz.G_post ) - Gz.t_ramp;  % Duration of pre-phasing gradient plateau [s]

tmp = 0 : sys.dt : Gz.t_post;  % Refocusing gradient gradient time vector [s]

Gz.G_post_ramp = ( 0 : -Gz.dG : Gz.G_post )';     % Refocusing gradient ramp waveform [T/m]
Gz.G_post_plat = ones( size(tmp') ) * Gz.G_post;  % Refocusing gradient plateau waveforms [T/m]

%% Composite the RF and Gradient waveforms and plot them
Gz.G = [0; Gz.G_ramp; Gz.G; flipud(Gz.G_ramp); Gz.G_post_ramp; Gz.G_post_plat; flipud(Gz.G_post_ramp); 0]; % Gradient amplitude [T/m]

% Timing corrections needed to center RF pulse
delay = zeros( length( Gz.G_ramp ), 1 );
pad = zeros( length( [Gz.G_post_ramp; Gz.G_post_plat; flipud(Gz.G_post_ramp); Gz.G_ramp] ), 1 );

RF.B1 = [0; delay; RF.B1; pad; 0]; % B1 amplitude [T]
Gx.G = zeros ( size ( Gz.G ) );   % Gradient amplitude [T/m]
Gy.G = zeros ( size ( Gz.G ) );   % Gradient amplitude [T/m]

%% Plot the final waveforms
Rad229_PSD_fig(1e6*RF.B1, Gx.G, Gy.G, Gz.G, sys.dt);

return
