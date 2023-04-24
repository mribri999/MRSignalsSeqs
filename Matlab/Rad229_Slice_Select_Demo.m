
%% Rad229_Slice_Select_Demo – Demonstrate designing both the RF pulse and gradient for slice selection.
% 
%  SYNTAX - [acq, sys, Gx, Gy, Gz, RF] = Rad229_Slice_Select_Demo(acq, RF);
%
%  INPUTS -  acq.dz = 10.0e-3;         % Slice thickness [m]
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

%% These questions can be used to further explore the code and concepts. Use the default settings as a starting point.
% RF Questions:
%     1.1) What is the highest possible flip angle without a B1 warning? 
% 
%     1.2) What is the highest possible TBW without a B1 warning?
%     
%     1.3) A 180' pulse exceeds B1,max. What adjustments can be made to enable a 180'?
%
%     1.4) Design the shortest possible TBW = 8 refocusing pulse that doesn't exceed B1,max.
%
%     1.5) [Advanced] Make a plot of B1,max and Gz,max as a function of flip angle for the shortest 
%          duration RF pulse. Use flip angles from 5 to 180 degrees. All pulses should be within ALL 
%          default hardware limits. What do you observe about the shape of the curve and the hardware
%          limits encountered?
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
% Slice Selection Questions
%     3.1) [Advanced] Compare the slice profiles for a 2ms TBW=4 and TBW=16 RF pulse. What differences 
%          do you observe?
%
%     3.2) [Advanced] Plot the slice profile with and without apodization.
%          What differences do you observe?
%
% Non-1H Imaging Questions (sys.gamma_bar=10.7084e6 [Hz/T]):
%     4.1) What happens if the default pulse is used to excite 13C? Why?
%
%     4.2) What else needs to change about the RF pulse to excite 13C?
%
%     4.3) [Advanced] Design an TBW=10 RF pulse to refocus 13C for a 10mm slice within hardware 
%          limits. If you were to build an MRI system for 13C imaging would you prefer an increase
%          in G_max or B1_max? Why?

function [acq, sys, Gx, Gy, Gz, RF] = Rad229_Slice_Select_Demo(acq, RF)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;
% sys.gamma_bar=10.7084e6;  % Gyromagnetic ratio for 11C

%% Define the RF pulse to excite a slice
if nargin == 0
  acq.dz = 10.0e-3;         % Slice thickness [m]

  RF.alpha = 90;            % RF flip angle [degrees]
  RF.TBW = 6;               % Time*bandwidth product [unitless] (relates to quality of RF pulse)
  RF.dur = 2.5e-3;          % RF pulse duration [seconds]
  RF.apod = 0.46;           % apod = 0 is non-apodized, apod = 0.5 is 'Hanning windowed', apod<0.5 is 'Hamming windowed'
end

% Parameters calculated from the INPUTS
RF.BW = RF.TBW ./ RF.dur;            % RF bandwidth [Hz]
RF.N = ceil( RF.dur / sys.dt );      % Number of points in the pulse
RF.t = linspace( -1, 1, RF.N )';     % Normalized time vector (duration is 1)

%% Design the RF pulse (waveform)
RF.apod_fxn = ( 1 - RF.apod ) + RF.apod * cos( pi * RF.t );  % RF pulse apodization function
RF.B1 = sinc( RF.t * RF.TBW / 2 ) .* RF.apod_fxn;            % Apodized (windowed) RF SINC pulse
RF.B1 = RF.B1 ./ sum( RF.B1 );                               % Normalize to unit area
RF.B1 = ( RF.alpha * pi/180 ) * RF.B1 / ( 2 * pi * sys.gamma_bar * sys.dt); % Scale for flip angle and real amplitude
  if max( RF.B1 ) > sys.B1max, warning('RF.B1 exceeds sys.B1max'); end

%% Design the gradient to excite a slice
Gz.G_amp = RF.BW / ( sys.gamma_bar * acq.dz ); % Gradient amplitude [T/m]
  if max(Gz.G_amp)>sys.G_max, warning('Gz.G_amp exceeds sys.G_max'); end
Gz.G_plat = Gz.G_amp * ones( size( RF.B1 ) );       % Gradient amplitude [T/m]

Gz.dG = sys.S_max * sys.dt;  % Slice-select gradient incremental gradient step [T/m]

Gz.G_ramp = ( 0 : Gz.dG : Gz.G_amp )'; % Readout gradient ramp [T/m]
%  if max(Gz.G_ramp) ~= Gz.G_amp, warning('Gz.G_ramp does not reach target value of Gz.G_amp'); end
Gz.t_ramp = sys.dt * numel(Gz.G_ramp);
% Gz.t_ramp = sys.G_max / sys.S_max;  % Ramp time [s]
%   if ~isinteger(Gz.t_ramp / sys.dt), warning('Gz.t_ramp is NOT an integer number of sys.dt!'); end

Gz.G = [0; Gz.G_ramp; Gz.G_plat; flipud(Gz.G_ramp)];
  
%% Design the post-excitation slice-select refocusing gradient (SSRG)
Gz.M0 = sum( Gz.G * sys.dt );  % The SSRG needs to have half this area [T*s/m]
Gz.G_ssrg_amp = -sys.G_max;    % The SSRG should use the maximum hardware to be fast [T/m] 
Gz.G_ssrg_ramp_time = sys.G_max ./ sys.S_max; % The SSRG should ramp as fast as possible

Gz.t_ssrg = ( ( Gz.M0 / 2 ) - 2 * (0.5) * Gz.G_ssrg_ramp_time * abs(Gz.G_ssrg_amp) ) / abs( Gz.G_ssrg_amp ); % Duration of pre-phasing gradient plateau [s]
  % The SSGR gradient needs half the moment, the plateau duration needs to account for moments from the two ramps too.

tmp = 0 : sys.dt : Gz.t_ssrg;  % Refocusing gradient gradient time vector [s]

Gz.G_ssrg_ramp = ( 0 : -Gz.dG : Gz.G_ssrg_amp )';     % Refocusing gradient ramp waveform [T/m]
Gz.G_ssrg_plat = ones( size(tmp') ) * Gz.G_ssrg_amp;  % Refocusing gradient plateau waveforms [T/m]

% tmp_M0 = 2 * sum( [ Gz.G_ssrg_ramp; Gz.G_ssrg_plat; flipud(Gz.G_ssrg_ramp) ] * sys.dt )

%% Composite the RF and Gradient waveforms and plot them
Gz.G = [0; Gz.G_ramp; Gz.G_plat; flipud(Gz.G_ramp); Gz.G_ssrg_ramp; Gz.G_ssrg_plat; flipud(Gz.G_ssrg_ramp); 0]; % Gradient amplitude [T/m]

% Timing corrections needed to center RF pulse
delay = zeros( length( Gz.G_ramp ), 1 ); % Delay RF start to end of gradient ramp
pad = zeros( length( [Gz.G_ssrg_ramp; Gz.G_ssrg_plat; Gz.G_ssrg_ramp; Gz.G_ramp] ), 1 );

RF.B1 = [0; delay; RF.B1; pad; 0]; % B1 amplitude [T]
Gx.G = zeros ( size ( Gz.G ) );    % Gradient amplitude [T/m]
Gy.G = zeros ( size ( Gz.G ) );    % Gradient amplitude [T/m]

%% Plot the final waveforms
% Rad229_PSD_fig(1e6*RF.B1, Gx.G, Gy.G, Gz.G, sys.dt);

return