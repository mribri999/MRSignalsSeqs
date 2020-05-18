%% Bval with Prisma gradients
% No limitation on PNS
params = struct;
params.mode = 'diff_bval';
params.gmax = 0.08;
params.smax = 200.0;
params.MMT = 2;
params.T_readout = 12.0;
params.T_90 = 3.0;
params.T_180 = 6.0;
params.dt = 100e-6;

G_min = get_min_TE_diff(100, 40.0, 80.0, params);
plot_waveform(G_min, params, 1, 1, 1)

%% Now with PNS = 1
% Note this is significantly slower because a convolution needs to be done
params = struct;
params.mode = 'diff_bval';
params.gmax = 0.08;
params.smax = 200.0;
params.MMT = 2;
params.T_readout = 12.0;
params.T_90 = 3.0;
params.T_180 = 6.0;
params.dt = 100e-6;
params.pns_thresh = 1.0;

G_min = get_min_TE_diff(100, 40.0, 80.0, params);
plot_waveform(G_min, params, 1, 1, 1)


%% Now just slew derated
% Not a huge difference, but PNS limited is better (and doesnt require
% figuring out the ideal smax)
params = struct;
params.mode = 'diff_bval';
params.gmax = 0.08;
params.smax = 50.0;
params.MMT = 2;
params.T_readout = 12.0;
params.T_90 = 3.0;
params.T_180 = 6.0;
params.dt = 200e-6;

G_min = get_min_TE_diff(100, 40.0, 80.0, params);
plot_waveform(G_min, params, 1, 1, 1)