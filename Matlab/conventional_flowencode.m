function [FE, DeltaM1] = conventional_flowencode(params)

    GAM = 2*pi*42.57e6;         % rad/(s*T)
    
    % M1 needed for flow encoding Venc
    DeltaM1 = pi/(GAM*params.VENC);  % T/m x s^2
    
    % Moments at end of slice select
    M0S = params.M0S(1,end);
    M1S = params.M1S(1,end);
   
    % Now solve for the next two lobes as presented in: 
    % Pelc NJ, Bernstein MA, Shimakawa A, Glover GH. Encoding strategies for three-direction phasecontrast MR imaging of flow. J Magn Reson Imaging 1991;1(4):405-413.
    % But we add an additional search over risetime (r) to allow for
    % waveforms with triangle lobes
    
    % Maximum riseitme [s]
    rmax = params.gmax/params.smax;
    
    % r = risetime in [s]
    for r = 1e-06:2e-06:rmax
        h = r*params.smax; % gradient strength [T/m]
        
        M02 = abs(-h*r+sqrt((h*r)^2+2*(h*r*M0S + M0S^2 + 2*h*(M1S+DeltaM1))))/2; % Area of the second lobe (see above paper for derivation) [T/m x s]
        M01 = M0S + M02; % Area of the first lobe (see above paper for derivation) [T/m x s]
        
        w1 = abs(M01)/h + r; % Duration of first lobe [s]
        w2 = abs(M02)/h + r; % Duration of second lobe [s]
        
        r_ = (ceil(r/params.dt)); % risetime in 'dt' units
        w1_ = (ceil(w1/params.dt)); % Duration of first lobe in 'dt' units
        w2_ = (ceil(w2/params.dt)); % Duration of second lobe in 'dt' units
        
        % We stop this loop if the flattop time is less than one 'dt' unit
        % (triangle lobe)
        if  (w1_-2*r_ <= 1) || (w2_-2*r_ <= 1)
            break
        end
        
    end
    
    % Build the flow compensated bipolar
    G = horzcat(linspace(0,-h,r_),linspace(-h,-h,w1_-2*r_),linspace(-h,0,r_),linspace(h/r_,h,r_-1),linspace(h,h,w2_-2*r_),linspace(h,0,r_));

    % Output waveform in [T/m]
    FE = horzcat(params.G_ss,G);

end
