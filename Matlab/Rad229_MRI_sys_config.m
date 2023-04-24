%% This function defines some default parameters and constants for an MRI system. 
% The MRI hardware is designed to be whimpy. This is helpful to better see the gradient 
% waveforms and to make the gradient waveform timing a little easier, thereby avoiding 
% some timing errors (maybe).
%
% DBE@STANFORD.EDU (March 2021) for Rad229

function sys = Rad229_MRI_sys_config
%% Default Parameters (don't change these!):
%  These are the default parameters that should be used to define and/or
%  reset the system. They are chosen specifically to help with making
%  certain phenomena more apparent.
%
% sys.B0 = 3.0; sys.B0_units = 'T';                            % Main (B0) field strength [T]
% sys.B1max = 25e-6; sys.B1max_units = 'T';                    % RF (B1) maximum field strength [T]
% sys.G_max = 10e-3; sys.G_max_unts = 'T/m';                   % Gradient maximum [T/m]
% sys.S_max = 100; sys.S_max_units = 'T/m/s';                  % Slewrate maximum [T/m/s]
% sys.dt = 10e-6; sys.dt_units = 's';                          % Waveform time steps (i.e. "raster time") [s]
% sys.gamma_bar=42.577478518e6; sys.gamma_bar_units = 'Hz/T';  % 1H gyromagnetic ratio [Hz/T]

%% Define the MRI system configuration
sys.B0 = 3.0; sys.B0_units = 'T';                            % Main (B0) field strength [T]
sys.B1max = 25e-6; sys.B1max_units = 'T';                    % RF (B1) maximum field strength [T]
sys.G_max = 10e-3; sys.G_max_unts = 'T/m';                   % Gradient maximum [T/m]
sys.S_max = 100; sys.S_max_units = 'T/m/s';                  % Slewrate maximum [T/m/s]
sys.dt = 10e-6; sys.dt_units = 's';                          % Waveform time steps (i.e. "raster time") [s]
sys.gamma_bar=42.577478518e6; sys.gamma_bar_units = 'Hz/T';  % 1H gyromagnetic ratio [Hz/T]
% sys.gamma_bar=10.7084e6; sys.gamma_bar_units = 'Hz/T';       % 13C gyromagnetic ratio [Hz/T]