function bval = get_bval(G, params)
%GET_BVAL Summary of this function goes here
%   Detailed explanation goes here

    if isfield(params, 'dt_out')
        dt = params.dt_out;
    else
        dt = params.dt;
    end
    
    T_READOUT = params.T_readout;
    
    TE = numel(G)*dt*1e3+T_READOUT;

    tINV = floor(TE/dt/1e3/2);

    GAMMA   = 42.58e3; 
    INV = ones(numel(G),1);   
    INV(tINV:end) = -1;

    Gt = cumsum(G'.*INV.*dt);
    bval = sum(Gt.*Gt);
    bval = bval * (GAMMA*2*pi)^2 * dt;
end

