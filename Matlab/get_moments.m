function moments = get_moments( G, T_READOUT, dt )
%GET_MOMENTS Summary of this function goes here
%   Detailed explanation goes here
    
    TE = numel(G)*dt*1e3+T_READOUT;

    tINV = floor(TE/dt/1e3/2);

    INV = ones(numel(G),1);   INV(tINV:end) = -1;
    
    Nm = 5;
    tvec = (0:numel(G)-1)*dt; % in sec
    % tMat for all moments
    tMat = zeros( Nm, numel(G) );
    for mm=1:Nm,
        tMat( mm, : ) = tvec.^(mm-1);
    end
    
    moments = abs(dt*tMat*(G'.*INV));
end

