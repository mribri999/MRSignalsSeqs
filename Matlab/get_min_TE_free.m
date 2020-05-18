function [G_out, T_out] = get_min_TE_free( params, T_hi )
%GET_MIN_TE Summary of this function goes here
%   Detailed explanation goes here

T_lo = 0.0;
T_range = T_hi-T_lo;

best_time = 999999.9;

if isfield(params, 'dt')
    dt = params.dt;
else
    dt = 1.0e-3/params.N0;
end

fprintf('Testing TE =');
while ((T_range*1e-3) > (dt/4.0))
    params.TE = T_lo + (T_range)/2.0;
    fprintf(' %.3f', params.TE);
    [G, lim_break, params] = gropt(params);
    if lim_break == 0
        T_hi = params.TE;
        if T_hi < best_time
            G_out = G;
            T_out = T_hi;
            best_time = T_hi;
        end
    else
        T_lo = params.TE;
    end
    T_range = T_hi-T_lo;
end

fprintf(' Final TE = %.3f ms\n', T_out);

end

