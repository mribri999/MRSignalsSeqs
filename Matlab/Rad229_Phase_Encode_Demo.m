%% Rad229_Phase_Encode_Demo – Demonstrate phase encode gradient design for imaging.
%
%  Recall, dy = FOVy / Ny . In general, the user can usually pick FOVy and Ny on the scanner.
%         dky = 1 / FOVy = gamma_bar * dG * t_PE . The k-space step size defines the FOV.
%
%  Note: See code for variable definitions and units
%  
% SYNTAX - [acq, sys, Gx, Gy, Gz, RF] = Rad229_Phase_Encode_Demo(acq)
%
% INPUTS -  acq.FOVy = 100e-3;  % Field-of-view along y-direction [m]
%           acq.Ny = 33;        % Number of pixels to discretize FOVy [#]
%           acq.y_pos = linspace(-acq.FOVy/2, acq.FOVy/2, acq.Ny); % Define the pixel locations [m]
%
% OUTPUTS – acq - A structure of the defined acquisition parameters.
%           sys - A structure of the defined MRI system parameters.
%           Gx - A structure with the Gx gradient parameters.
%           Gy - A structure with the Gy gradient parameters.
%           Gz - A structure with the Gz gradient parameters.
%           RF - A structure with the RF gradient parameters.   
%
% DBE@STANFORD.EDU (March 2021) for Rad229

%% These questions can be used to further explore the code and concepts.
%  Phase Encode Questions:
%     1) How many cycles are there across the FOV when using the maximum phase encoding gradient?
%        Does this satisfy Nyquist?
%
%     2) How many cycles per pixel? Decrease the resolution and answer the same questions. Discuss.
%     
%     3) What happens to the duration of the phase encoding gradients for higher resolution imaging?
%
%     4) [Advanced] Why are the phase encoding gradient timing warnings important? If there is a 
%        warning, then what would be the artifact? What causes the possible timing error?
%     
%     5) [Advanced] Incorporate a FOV shift in the PE direction. Hint: See Bernstein p. 264

function [acq, sys, Gx, Gy, Gz, RF] = Rad229_Phase_Encode_Demo(acq)

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Define Phase Encoding Variables - Most scanners are configured for the user to select the FOV and 
%  matrix size for encoding, which then defines the spatial resolution.
if nargin == 0
  acq.FOVy = 100e-3;  % Field-of-view along y-direction [m]
  acq.Ny = 33;        % Number of pixels to discretize FOVy [#]
%   acq.Ny = 65;        % Number of pixels to discretize FOVy [#]
  acq.y_pos = linspace(-acq.FOVy/2, acq.FOVy/2, acq.Ny); % Define the pixel locations [m]
end
acq.dy = acq.FOVy / acq.Ny;     % Pixel dimension along y-direction [m]

%% Design the phase encode gradients
Gy.dky = 1 / acq.FOVy;                    % Spatial frequency increment [cycle/m]
Gy.ky_max = ( (acq.Ny-1) / 2 ) * Gy.dky;  % Maximum phase encode spatial frequency [cycle/m]

% Phase encode gradient timing [ramp, then plateau]
Gy.t_ramp = sys.G_max / sys.S_max;  % Ramp time [s]
%   if ~isinteger(Gy.t_ramp / sys.dt), warning('Gy.t_ramp integer number of sys.dt!'); end

% Recall, k = gamma_bar * G * t
Gy.t_plat = ( Gy.ky_max / (sys.gamma_bar * sys.G_max) ) - Gy.t_ramp;  % Plateau time [s]
%   if ~isinteger(Gy.t_plat / sys.dt), warning('Gy.t_plat integer number of sys.dt!'); end
  if Gy.t_plat<0, error('Gy.t_plat is NEGATIVE!'); end

% Phase encode gradient waveform design
Gy.dG = Gy.dky / (sys.gamma_bar * (Gy.t_ramp + Gy.t_plat) ); % Gradient amplitude increment [T/m]
Gy.G_steps = -sys.G_max : Gy.dG : sys.G_max;  % Phase encode gradient amplitudes [T/m]
tmp = 0 : sys.dt : Gy.t_plat;  % Phase encode gradient time vector [s]
Gy.G_plat = ones( size(tmp') ) * Gy.G_steps;  % Phase encode plateau waveforms [T/m]

% Design each phase encoding step
for n = 1 : acq.Ny
  m = Gy.t_ramp / sys.dt;  % Number of raster steps along phase encode ramp
  grad_step = abs( Gy.G_steps(n) / m ) + eps; % Phase encode gradient increment (unique for each PE step)
  Gy.G_ramp(:,n) = sign(Gy.G_steps(n)) * ( 0 : grad_step : abs( Gy.G_steps(n) ) );
    % Design the precise gradient raster amplitudes (bookkeeping to get sign right)
end

% Composite the first ramp, plateau, and last ramp
Gy.G = [Gy.G_ramp; Gy.G_plat; flipud(Gy.G_ramp) ];

%% Plot the waveforms and visualize the encoding per FOV and per pixel
% Define a subset of PE steps to examine (cubic sampling looks good...)
m = (acq.Ny-1) / 2;
% PE_IND = unique(round(m * ( linspace(-m, m, 7) ).^3 / (m^3) + (m+1) ) ); % Overly complicated...
PE_IND = 1 : 4 : acq.Ny; % Used during plotting to show a subset of PE steps

% Define the RF and gradient waveforms
RF.B1 = zeros( size(Gy.G , 1) , 1 );
Gx.G  = zeros( size(Gy.G , 1) , 1 );
Gy.G  = Gy.G;
Gz.G  = zeros( size(Gy.G , 1) , 1 );

% Make a nice figure
Rad229_PSD_fig(1e6*RF.B1, Gz.G, Gy.G(:,PE_IND), Gx.G, sys.dt);

%% What do the phase encode gradients do to the phase of the transverse magnetization in each pixel?
pix = interp(acq.y_pos,10); % Upsample the positions used to sample the encoding (makes visualization a bit better)

figure; hold on; 
IND = 1;
for k = PE_IND
  phase = (2*pi) * sys.gamma_bar * sum( Gy.G(:,k) * sys.dt ) * pix; % Bulk magnetization phase depends on gradient area and position
  subplot(length(PE_IND),1,IND); hold on;
  q = plot([acq.y_pos(1:10:end); acq.y_pos(1:10:end)],[-1 1],'k'); set(q,'LineWidth',2); % Plot some reference position lines
  p = plot(pix,sin(phase)); set(p,'LineWidth',3); axis tight;
  if IND == 1, title('Phase Encode Gradient Induced Phase Across FOV'); end
  if IND == length(PE_IND), xlabel('Field Of View [m]'); end
  IND=IND+1;
end

return