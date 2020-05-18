%% Design a slice select gradient
Nm = 3; % Number of moments to calculate
dt = 10e-6;

% Hand tuned to get M0ss and M1ss roughly = 1.4, .133
g_ss = .0070;
N_ss = round(1.375/(g_ss*dt)*1e-6);
G_ss = ones(1, N_ss) .* g_ss;
t_ss = numel(G_ss) * dt;

tvec = (0:numel(G_ss)-1)*dt*1e3; % in msec
% tMat for all moments
tMat = zeros( Nm, numel(G_ss) );
for mm=1:Nm,
    tMat( mm, : ) = tvec.^(mm-1);
end
moments_ss = (1e3*1e3*dt*tMat*(G_ss'));

M0ss = moments_ss(1);
M1ss = moments_ss(2);
%% Compute flow comped refocuser
% note that t_ss has been added to the third spot in moment params
% this tells the optimizer that the first point in the waveform is actually
% at t = t_tss instead of 0
params = struct;
params.mode = 'free';
params.gmax = 0.05;
params.smax = 100.0;
params.moment_params = [];
params.moment_params(:,end+1) = [0, 0, t_ss, -1, -1, -M0ss, 1.0e-4];
params.moment_params(:,end+1) = [0, 1, t_ss, -1, -1, -M1ss, 1.0e-4];
params.dt = dt;

[G_min, T_min] = get_min_TE_free(params, 1.0);

%%
G = [G_ss G_min];
figure()
plot(G)

tvec = (0:numel(G)-1)*dt*1e3; 
tMat = zeros( Nm, numel(G) );
for mm=1:Nm,
    tMat( mm, : ) = tvec.^(mm-1);
end

moments = (1e3*1e3*dt*tMat*(G'));
fprintf('Final moments = %.03f  %.03f\n', moments(1), moments(2));