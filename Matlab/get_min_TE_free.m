function [G_out, T_out] = get_min_TE_free( params, T_hi )
%GET_MIN_TE Finds the shortest possible timing for the given GrOpt parameters
%
% SYNTAX  - [G_out, T_out] = get_min_TE_free( params, T_hi )
%
% INPUTS  - params: the parameters that go to the GrOpt(params) call [struct]
%         - T_hi: maxmimum gradient waveform duration to search [float, ms]
%
% OUTPUTS - G_out: Output gradient waveform [Nx1 float array, T/m]
%           T_out: Duration of G_out (this is simply numel(G) * dt) [float, ms]
%
% The initial temporal search duration spans [0, T_hi] in ms.
% (Currently this just equals T_hi, but there may be situations where it
% could be better to have nonzero T_lo)
%
% MLoecher@STANFORD.EDU (March 2020) original code.
% DBE@STANFORD.EDU (March 2021 updates)

T_lo = 0.0;
T_range = T_hi-T_lo;

best_time = 999999.9; % We store the optimal waveform duration here, in ms

if isfield(params, 'dt')
    dt = params.dt;  % Raster time in s
else
    dt = 1.0e-3/params.N0;  % Raster time in s
end

fprintf('Testing TE =');

% We do a bisecting line search, and stop when the section is smaller than
% dt/4 (therefore no more change in waveform size N)
while ((T_range*1e-3) > (dt/4.0))
    params.TE = T_lo + (T_range)/2.0;  % TE is the waveform duration in ms, calculated by bisecting (T_lo,T_hi)
    fprintf(' %.3f', params.TE);
    [G, lim_break, params] = gropt(params); % Run GrOpt on params
    if lim_break == 0  % lim_break == 0 means all constraints are met, so the waveform is good
        T_hi = params.TE;  % We know the current TE is good, so it becomes the upper bound of our search
        if T_hi < best_time  % If this is the best time (it should be), update the best_time and set best outputs
            G_out = G; 
            T_out = T_hi;
            best_time = T_hi;
        end
    else 
        % If the waveform does not pass constraints, it is too short, so set the lower bound of our 
        % search range to this time (because everything shorter should also fail)
        T_lo = params.TE;
    end
    T_range = T_hi-T_lo;  % This is just used to calculate the next bisection
end

fprintf(' Final TE = %.3f ms\n', T_out);

end