function [TE,G,b] = design_symmetric_gradients(bvalue_target,T_ECHO,T_90,G_Max,MMT)
% Returns the TE for symmetric DWI waveforms with a specified b-value and
% sequence timing parameters. The waveforms used are: MONOPOLAR, BIPOLAR
% and MODIFIED BIPOLAR (Stoeck CT, von Deuster C, Genet M, Atkinson D, Kozerke
%                       S. Second-order motion-compensated spin echo diffusion
%                       tensor imaging of the human heart. MRM. 2015.)
%
% INPUTS:  G_Max        - Max gradient amplitude [T/m]
%          bvalue_target- Desired b-value [s/mm2]
%          T_ECHO       - EPI time to Echo [ms]
%          T_90         - Start time of diffusion encoding [ms]
%          MMT          - Desired waveform moments
%                         - 0 - M0= 0      - MONO
%                         - 1 - M0M1 = 0   - BIPOLAR
%                         - 2 - M0M1M2 = 0 - MODIFIED BIPOLAR
%
% OUTPUTS: TE    -  TE of resultant waveform [ms]
%          G     -  Diffusion encoding waveform [T/m]
%          b     -  b-value of encoding waveform [s/mm2]
%
% Magnetic Resonance Research Labs (http://mrrl.ucla.edu)
% Department of Radiological Sciences
% University of California, Los Angeles
% Eric Aliotta (EAliotta@mednet.ucla.edu)
% Holden Wu (HWu@mednet.ucla.edu)
% Daniel Ennis (DEnnis@mednet.ucla.edu)
% December 16, 2015


epsilon = 1.5;   % gradient ramp time 
RFgap     = 4.3; % 180 pulse duration
epsilon = floor(epsilon*10)/10;

% define monopolar waveform
if MMT == 0
  gap = RFgap;
  N = 4*epsilon + gap + 2*T_ECHO+T_90; % starting point for total duration
  T = 0.1; % scale time in ms
  b = 0;
  
  % update waveform until the b-value is large enough
  while(b<bvalue_target*0.995)
    N = N+T;
    
    time = N;
    
    lambda = (time - 4*epsilon - gap - 2*T_ECHO)/2;
    lambda = round(lambda/T);
    
    grad = trapTransform([lambda,lambda],G_Max,floor(epsilon/T),1,floor((T_ECHO-T_90+gap)/T),1);
    
    n = length(grad);
    
    C=tril(ones(n));
    C2 = C'*C;
    GAMMA = 42580;
    
    INV = ones(n,1);   INV(floor((n+T_ECHO)/2):end) = -1;
    
    Ts = T*(1e-3);
    
    b = (GAMMA*2*pi)^2*(grad.*INV*Ts)'*(C2*(grad.*INV*Ts))*Ts;
    
    tINV = ceil(lambda + floor((T_ECHO-T_90+gap)/T) + 2*epsilon/T - 0.5*gap/T);
    TEh1 = T_ECHO/T + length(grad(tINV:end));
    TEh2 = tINV;
    
    TE = 2*max(TEh1,TEh2)*T;
    G = grad;
  end
end

% define bipolar waveform (M1=0)
if MMT == 1
  
  L = 1; % starting point
  T = 0.1; % scale time in ms
  
  b = 0;
  % update waveform until the b-value is large enough
  while(b<bvalue_target*0.995)
    
    L = L+T;
    
    % duration of each bipolar lobe
    lambda  = L;         
    LAMBDA  = lambda;    
    
    LAMBDA = round(LAMBDA/T);
    lambda = round(lambda/T);

    % gap b/w gradients is just the RF duration 
    gap = RFgap;
    
    % take trapezoid durations and create G(t) vector
    grad = trapTransform([lambda,-LAMBDA,LAMBDA,-lambda],G_Max,round(epsilon/T),1,round(gap/T),2);
    
    % legnth of waveform
    n = length(grad);
    
    % vector for b-value integration
    C=tril(ones(n));
    C2 = C'*C;
    
    % Gyromagnetic ratio
    GAMMA = 42580;
    
    % refocusing pulse time
    tINV = floor(n/2);
    
    % vector to invert magnetization (+1 before 180, -1 after 180)
    INV = ones(n,1);
    INV(floor(tINV):end) = -1;
    
    % time increment in seconds
    Ts = T*(1e-3);
    
    % test b-value
    b = (GAMMA*2*pi)^2*(grad.*INV*Ts)'*(C2*(grad.*INV*Ts))*Ts;
    
    % pre 180 contribution to TE
    TEh1 = 0.5*RFgap/T + lambda + LAMBDA + 4*epsilon/T + T_ECHO/T;
    % post 180 contribution to TE
    TEh2 = 0.5*RFgap/T + lambda + LAMBDA + 4*epsilon/T + T_90/T;
    
    % Final TE
    TE = 2*max(TEh1,TEh2)*T + 2 + 2; %additional 4ms for spoilers.
    
    % final gradient waveform
    G = grad;
  end
end

% define modified bipolar (M1=M2 = 0) waveform
if MMT == 2
  L = 1; % starting point
  T = 0.1; % scale in ms
  
  b = 0;
  
  % update waveform until the b-value is large enough
  while(b<bvalue_target*0.995)
    
    L = L+T;
    
    % first trap duration
    lambda = L;                     lambda = round(lambda/T); 
    % second trap duration
    LAMBDA  = 2*lambda + epsilon;   LAMBDA = round(LAMBDA/T);
    
    % time between first and second sets of gradients
    gap = 2*epsilon + lambda;
    
    % take trapezoid durations and create G(t) vector
    grad = trapTransform([lambda,-LAMBDA,-LAMBDA,lambda],G_Max,round(epsilon/T),1,round(gap/T),2);
    
    % legnth of waveform
    n = length(grad);
    
    % vector for b-value integration
    C=tril(ones(n));
    C2 = C'*C;
    
    % Gyromagnetic ratio
    GAMMA = 42580;
    
    % refocusing pulse time
    tINV = n/2 + round(gap/T) - round(RFgap/T);
    
    % vector to invert magnetization (+1 before 180, -1 after 180)
    INV = ones(n,1);
    INV(floor(tINV):end) = -1;
    
    % time increment in seconds
    Ts = T*(1e-3);
    
    % test b-value
    b = (GAMMA*2*pi)^2*(grad.*INV*Ts)'*(C2*(grad.*INV*Ts))*Ts;
    
    % pre 180 contribution to TE
    TEh1 = 0.5*RFgap/T + lambda + LAMBDA + 4*epsilon/T + T_ECHO/T;
    % post 180 contribution to TE
    TEh2 = -0.5*RFgap/T + lambda + LAMBDA + 4*epsilon/T + T_90/T + gap/T;
    
    % final TE
    TE = 2*max(TEh1,TEh2)*T;
    
    % final gradient waveform
    G = grad;
  end
end


end

function g = trapTransform(f,G0,SR,tFact,gap,gpos)
% gradient waveform from trapezoidal reduction
% define waveform in terms of gradient duration and sign assuming 
% G = Gmax 
%
% input: f  -- row of numbers indicating the duration of each gradient lobe
%              in ms. Must correspond to an integer number of timepoints
%        G0 -- Gmax. All lobes assumed to be at Gmax
%        SR -- Slew duration (normalized to unit size) (default- 1)
%        tFact  -- Temporal resolution subsampling (default- 1)
%        gap-- Gap duration an RF pulse [units] (default 0)
%        gpos- Position of gap (list the number of f entry to put the gap
%              AFTER (default floor(num entries/2)
% 
% output: g -- fully represented gradient waveform

if nargin<2
    G0 = 0.074;
end

if nargin<3
    SR = 1;
end

if nargin<4
    tFact = 1;
end

if nargin<5
    gap = 0;
end

if nargin<6
    gpos = floor(length(f)/2);
end

if tFact == 1e-3
    tFact = 1;
    fprintf('Assuming T = 1ms, subsample of 1 used!! \n');
end

if min(abs(f)) < 1
    fprintf('ERROR - Need to allow time for slewing!!!\n');
    return;
end

%g = G0*ones( (sum(abs(f))+gap )*tFact,1);
g = G0*ones( (sum(abs(f)) + gap + 2*numel(f)*SR - (numel(f)-1) )*tFact,1);

count = 1;

for j=1:length(f)
    PLAT = abs(f(j));
    if j == gpos
        tnow = count;
        % ramp up
        g(tnow:tnow+SR-1) = g(tnow:tnow+SR-1).*(0:1/SR:1-1/SR)'*(f(j)/PLAT);
        %g(count+1:count+1+SR) = g(count+1:count+1+SR).*(0:1/SR:1)'*(f(j)/PLAT);
        tnow = tnow + SR;
        % plateau
        g(tnow:tnow+PLAT*tFact-1) = g(tnow:tnow+PLAT*tFact-1)*f(j)/PLAT;
        %g(count+2+SR:count+PLAT*tFact-SR+1) = g(count+2+SR:count+PLAT*tFact-SR+1)*f(j)/PLAT;
        tnow = tnow + PLAT*tFact;
        % ramp down
        g(tnow:tnow+SR-1) = g(tnow:tnow+SR-1).*(1-(1/SR:1/SR:1))'*(f(j)/PLAT);
       
        count = tnow + SR-1;
        
        g(count+1:count+gap*tFact) = g(count+1:count+gap*tFact)*0;
        count = count + gap*tFact;
    else
        tnow = count;
        % ramp up
        g(tnow:tnow+SR-1) = g(tnow:tnow+SR-1).*(0:1/SR:1-1/SR)'*(f(j)/PLAT);
        %g(count+1:count+1+SR) = g(count+1:count+1+SR).*(0:1/SR:1)'*(f(j)/PLAT);
        tnow = tnow + SR;
        % plateau
        g(tnow:tnow+PLAT*tFact-1) = g(tnow:tnow+PLAT*tFact-1)*f(j)/PLAT;
        %g(count+2+SR:count+PLAT*tFact-SR+1) = g(count+2+SR:count+PLAT*tFact-SR+1)*f(j)/PLAT;
        tnow = tnow + PLAT*tFact;
        % ramp down
        g(tnow:tnow+SR-1) = g(tnow:tnow+SR-1).*(1-(1/SR:1/SR:1))'*(f(j)/PLAT);
       count = tnow + SR-1;
        % ramp up
%         g(count+1:count+1+SR) = g(count+1:count+1+SR).*(0:1/SR:1)'*(f(j)/PLAT);
%         % plateau
%         g(count+1+SR:count+PLAT*tFact-SR) = g(count+1+SR:count+PLAT*tFact-SR)*f(j)/PLAT;
%         % ramp down
%         g(count+PLAT*tFact-SR:count+PLAT*tFact) = g(count+PLAT*tFact-SR:count+PLAT*tFact).*(1-(0:1/SR:1))'*(f(j)/PLAT);
%         
        %count = count + PLAT*tFact;
    end
end

end
