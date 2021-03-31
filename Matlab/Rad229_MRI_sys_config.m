%% This function defines some default parameters and constants. The MRI
% hardware is designed to be whimpy. This is helpful to better see the
% gradient waveforms and to make the gradient waveform timing a little
% easier, thereby avoiding timing errors.
%
% DBE@STANFORD.EDU (March 2021) for Rad229

function sys = Rad229_MRI_sys_config
%% Default Parameters:
%     sys.B0 = 3.0;             % Main (B0) field strength [T]
%     sys.B1max = 20e-6;        % RF (B1) maximum field strength [T]
%     sys.G_max = 10e-3;        % Gradient maximum [T/m]
%     sys.S_max = 100;           % Slewrate maximum [T/m/s]
%     sys.dt = 10e-6;           % Waveform time steps [s]
%     sys.gamma_bar=42.57e6;    % 1H gyromagnetic ratio [Hz/T]

%% Define the MRI system configuration
sys.B0 = 3.0;                  % Main (B0) field strength [T]
sys.B1max = 20e-6;             % RF (B1) maximum field strength [T]
sys.G_max = 10e-3;             % Gradient maximum [T/m]
sys.S_max = 100;               % Slewrate maximum [T/m/s]
sys.dt = 10e-6;                % Waveform time steps [s]
sys.gamma_bar=42.577478518e6;  % 1H gyromagnetic ratio [Hz/T]
% sys.gamma_bar=10.7084e6;     % 13C gyromagnetic ratio [Hz/T]