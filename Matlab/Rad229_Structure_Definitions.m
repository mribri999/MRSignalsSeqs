%% This is a README file that contains some information about how STRUCTURES are used to contain 
%  variables for the Rad229_* functions. The following base structures are
%  defined:
%
% sys.* for MRI system hardware, field, and imperfection specifications.
% 
%
% NOTE - This code is not meant to be run per se, but it does provide simple examples.
%
% DBE@STANFORD.EDU (May 2021) for Rad229

%% Define the MRI system configuration
sys.B0 = 3.0; sys.B0_units = 'T';                            % Main (B0) field strength [T]
sys.B1max = 20e-6; sys.B1max_units = 'T';                    % RF (B1) maximum field strength [T]
sys.G_max = 10e-3; sys.G_max_units = 'T/m';                   % Gradient maximum [T/m]
sys.S_max = 100; sys.S_max_units = 'T/m/s';                  % Slewrate maximum [T/m/s]
sys.dt = 10e-6; sys.dt_units = 's';                          % Waveform time steps (i.e. "raster time") [s]
sys.gamma_bar=42.577478518e6; sys.gamma_bar_units = 'Hz/T';  % 1H gyromagnetic ratio [Hz/T]

sys.t_ramp_max = sys.G_max ./ sys.S_max; % Maximum slew duration for any ramp [s]
sys.dG_max = sys.S_max * sys.dt;  % Maximum gradient increment [T/m]

% sys.off_x = linspace( 0 , 5 , acq.Nx ) .^ 2;  % Arbitrary off-resonance frequency [Hz]
% sys.off_y = linspace( 0 , 5 , acq.Nx ) .^ 2;  % Arbitrary off-resonance frequency [Hz]
% sys.off_xy = sys.off_y' * sys.off_x;    % Arbitrary off-resonance frequency across the FOV

% %% Acquisition parameters
% acq.Nx = 128;                 % Number of x-axis encoding points [#]
% acq.Nx = 128;                 % Number of y-axis encoding points [#]
% acq.ESP = 500e-6;             % Echo spacing [s]
% 
% %% RF Pulse parameters
% RF.alpha = 90;                       % RF flip angle [degrees]
% RF.TBW = 6;                          % Time*bandwidth product [unitless] (relates to quality of RF pulse)
% RF.dur = 2.0e-3;                     % RF pulse duration [seconds]
% RF.apod = 0.46;                      % apod = 0 is non-apodized, apod = 0.5 is 'Hanning windowed', apod<0.5 is 'Hamming windowed'
% 
% % Calculated parameters
% RF.BW = RF.TBW ./ RF.dur;            % RF bandwidth [Hz]
% RF.N = ceil( RF.dur / sys.dt );      % Number of points in the pulse
% RF.t = linspace( -1, 1, RF.N )';     % Normalized time vector (duration is 1)
% RF.apod_fxn = ( 1 - RF.apod ) + RF.apod * cos( pi * RF.t );  % RF pulse apodization function
% RF.B1 = sinc( RF.t * RF.TBW / 2 ) .* RF.apod_fxn;            % Apodized (windowed) RF SINC pulse
% RF.B1 = RF.B1 ./ sum( RF.B1 );                               % Normalize to unit area
% RF.B1 = ( RF.alpha * pi/180 ) * RF.B1 / ( 2 * pi * sys.gamma_bar * sys.dt); % Scale for flip angle and real amplitude
% 
% %% Object parameters
% % Chemical shift parameters
% obj.fat_ppm = -3.3e-6;                    % Chemical shift of fat [ppm];
% % obj.silicone_ppm = +4.5e-6;               % Chemical shift of silicone implant [ppm]; https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.22872
% obj.dB0_fat = obj.fat_ppm * sys.B0;       % B0 shift [T]
% obj.f_fat = sys.gamma_bar * obj.dB0_fat;  % Frequency offset [Hz]
% obj.num = 3; % Identify the phantom sub-object for fat (default phantom has 11 objects)

%% Gradient Objects - For conventional Cartesian imaging most gradient events can be represented as
%  a series of ramps and plateaus for each gradient axis. Hence we define:

% Define each gradient event
Gx(1).G_start = 0;
Gx(1).G_end = sys.G_max / pi; % Arbitrary gradient amplitude [T/m]
Gx(1).dur = sys.t_ramp_max;
Gx(1).name = 'ramp up';

Gx(2).G_start = Gx(1).G_end;
Gx(2).G_end = Gx(1).G_end;
Gx(2).dur = 217 * sys.dt; % Arbitrary gradient plateau duration [s]
Gx(2).name = 'plateau';

Gx(3).G_start = Gx(1).G_end;
Gx(3).G_end = 0; % Arbitrary gradient amplitude [T/m]
Gx(3).dur = sys.t_ramp_max;
Gx(4).name = 'ramp down';

% Generate each gradient waveform segment
for n = 1 : numel( Gx )
  Gx(n).N_pts = ( Gx(n).dur ./ sys.dt ) ;
  Gx(n).G = linspace( Gx(n).G_start , Gx(n).G_end , Gx(n).N_pts );
end

% I don't really like adding a "final" waveform that is the concat of the
% elemental waveforms...
% Combine the gradient waveform segments
Gx(end+1).G = cat(2,Gx.G);  % The gradient waveform for the FIRST x-gradient event