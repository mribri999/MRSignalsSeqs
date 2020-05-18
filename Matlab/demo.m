%% Simple diffusion solver
params.mode = 'diff_bval';
params.gmax = 0.04;
params.smax = 200.0;
params.MMT = 1;
params.TE = 60.0;
params.T_readout = 12.0;
params.T_90 = 3.0;
params.T_180 = 6.0;
params.dt = 500e-6;
params.dt_out = 10e-6;

[G, ~] = gropt(params);

plot_waveform(G, params)

%% Diffusion min_TE finder
G_min = get_min_TE_diff(250, 30.0, 200.0, params);

plot_waveform(G_min, params)

%% Bipolar generator

params.mode = 'free';
params.gmax = 0.04;
params.smax = 200.0;
% The structure of moment params entries is:
% [axis dir, moment order, 0 reference time, 
%        start time, end time, desired moment, tolerance]
params.moment_params = [];
params.moment_params(:,end+1) = [0, 0, 0, -1, -1, 0, 1.0e-3];
params.moment_params(:,end+1) = [0, 1, 0, -1, -1, 11.74, 1.0e-3];
params.TE = 1.32;
params.dt = 10e-6;

[G, ~] = gropt(params);

plot(G)

%% TE finder for free mode
% The only needed input parameter is the max search range

[G_min, T_min] = get_min_TE_free(params, 3.0);
plot(G_min)