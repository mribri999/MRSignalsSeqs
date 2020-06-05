function [s,P] = epg_X_rfspoil(flipangle,rf_phase,T1,T2,TR,f,ka,WT)
% Simulate RF spoiled sequence, 
% flipangle: flipangle in radians
% rf_phase: the vector of phases for the sequence of RF pulses
% T1: 1x2 vector with the T1 constant of both comparment
% T2: 1x2 vector or scalar with the T2 constants
% ka: forward exchange rate
% f: fraction of second compartment
% WT: saturation of the bound pool componen will be exp(-WT), ignore for
% full case

if exist('WT','var')
    P = epg_X_m0(f,'MT');
    caso = 'MT';
else
    P = epg_X_m0(f,'full');
    caso = 'full';
end

P = [P zeros(size(P,1),250)];
N = length(rf_phase);
for n=1:N
    P = epg_X_relax(P,TR,T1,T2,ka,0,f);
    P = epg_X_grad(P,1);  
  if strcmp(caso,'MT')
    P = epg_X_rf(P,flipangle,rf_phase(n),WT);	% RF excitation
        s(1,n) =  P(1,1);
  else
      P = epg_X_rf(P,flipangle,rf_phase(n));	% RF excitation
        s(1:2,n) = P([1 4],1);  	% Signal is F0 state.
  end
end
