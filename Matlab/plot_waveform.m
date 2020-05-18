function plot_waveform( G, params, plot_moments, plot_slew, plot_pns)
%PLOT_WAVEFORM Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    plot_moments = 1;
end
if nargin < 5
    plot_slew = 0;
end
if nargin < 6
    plot_pns = 0;
end

if isfield(params, 'dt_out')
    dt = params.dt_out;
else
    dt = params.dt;
end

if isfield(params, 'T_readout')
    T_readout = params.T_readout;
else
    T_readout = 0.0;
end


n_plots = 1;
n_plots = n_plots + plot_moments;
n_plots = n_plots + plot_slew;
n_plots = n_plots + plot_pns;

figure();

bval = get_bval(G, params);
TE2 = numel(G)*dt*1e3+T_readout;


% --------------
% Plot waveform
ii = 1;
subplot(n_plots,1,ii);
x = [0:numel(G)-1] .* dt * 1e3;
y = G * 1000;
plot(x,y,'LineWidth', 2);
ylim([-1.1 * max(abs(y)), 1.1 * max(abs(y))])
ylabel('G [mT/m]')
xlabel('t [ms]')
title(['bval = ' num2str(bval) '   TE = ' num2str(TE2)])
ii = ii + 1;


% --------------
% Plot moments
if plot_moments > 0
    tINV = floor(TE2/dt/1e3/2);

    INV = ones(numel(G),1);   INV(tINV:end) = -1;
    INV = INV';

    Nm = 5;
    tvec = (0:numel(G)-1)*dt; % in sec
    tMat = zeros( Nm, numel(G) );
    for mm=1:Nm,
        tMat( mm, : ) = tvec.^(mm-1);
    end

    moments = dt*tMat.*repmat(G.*INV, size(tMat,1), 1);
    moments = cumsum(moments, 2);

    subplot(n_plots,1,ii);
    hold on
    for i = 1:3
        plot(x, moments(i,:)/max(abs(moments(i,:))), 'LineWidth', 2);
    end

    ylabel('Moments [A.U]')
    xlabel('t [ms]')
    legend('M0', 'M1', 'M2')
    hline = refline([0 0]);
    hline.Color = 'k';
    ii = ii + 1;
end

% --------------
% Plot slew
if plot_slew > 0
    subplot(n_plots,1,ii);
    y = diff(G) / dt;
    y = [y 0];
    plot(x,y,'LineWidth', 2);
    ylabel('Slewrate [T/m/s]')
    xlabel('t [ms]')
    ii = ii + 1;
end

% --------------
% Plot slew
if plot_pns > 0
    subplot(n_plots,1,ii);
    y = abs(get_stim(G, dt)');
    plot(x,y,'LineWidth', 2);
    ylabel('PNS [A.U.]');
    xlabel('t [ms]');
    hline = refline([0 1.0]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    ii = ii + 1;
end

end

