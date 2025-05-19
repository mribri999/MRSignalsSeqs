function res = get_min_TE_diff( target_bval, min_TE, max_TE, params )
%GET_MIN_TE Summary of this function goes here
%   Detailed explanation goes here

if isfield(params, 'dt_out')
    dt = params.dt_out;
else
    dt = params.dt;
end

a = min_TE;
b = max_TE;
h = b - a;
tol = params.dt*1e3;

invphi = (sqrt(5) - 1) / 2;                                                                                                      
invphi2 = (3 - sqrt(5)) / 2;

n = ceil(log(tol/h)/log(invphi));

c = a + invphi2 * h;
d = a + invphi * h;

% disp([a,b,c,d]);

params.TE = c;
Gc = gropt(params);
yc = abs(get_bval(Gc, params) - target_bval);

params.TE = d;
Gd = gropt(params);
yd = abs(get_bval(Gd, params) - target_bval);

fprintf('Finding TE (%d iterations): ', n);
for k = 1:n
    fprintf('%d ', k);
    if (yc < yd)
        b = d;
        d = c;
        yd = yc;
        h = invphi*h;
        c = a + invphi2 * h;
        params.TE = c;
        Gc = gropt(params);
        yc = abs(get_bval(Gc, params) - target_bval);
    else
        a = c;
        c = d;
        yc = yd;
        h = invphi*h;
        d = a + invphi * h;
        params.TE = d;
        Gd = gropt(params);
        yd = abs(get_bval(Gd, params) - target_bval);
    end
end

% bvalc = get_bval(Gc, T_readout, dt);
% bvald = get_bval(Gd, T_readout, dt);
% if (bvalc < bvald)
%     res = Gd;
% else
%     res = Gc;
% end

if (yc > yd)
    res = Gc;
else
    res = Gd;
end

bval_out = get_bval(res, params);
TE_end = numel(res)*dt*1e3+params.T_readout;
if ( get_bval(res, params) < 0.90 * target_bval )
    fprintf(' TE = %f  bval = %f not big enough, restarting search\n', TE_end, bval_out );
    res = get_min_TE_diff( target_bval, 0.9*TE_end, 2*TE_end, params);
end

fprintf(' Done  TE = %f\n', numel(res)*dt*1e3+params.T_readout);

end