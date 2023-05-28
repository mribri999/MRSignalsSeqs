function [s,P] = epg_X_CMPG(flipangle,f,T1,T2,esp,ka,deltab)
%Simulate EPG-X CMPG sequence. Only for the T2 case, not the MT case
% flipangle: flipangle in radians
% T1: 1x2 vector with the T1 constant of both comparment
% T2: 1x2 vector or scalar with the T2 constants
% ka: forward exchange rate
% f: fraction of second compartment

etl = length(flipangle);
P = epg_X_m0(f,'full');
P = [P zeros(6,2*etl)];

% -- 90 excitation
P = epg_X_rf(P,pi/2,pi/2);	% Do 90 tip.
s = zeros(3,etl);		% Allocate signal vector to store.

for ech=1:etl
  P = epg_X_relax(P,esp/2,T1,T2,ka,deltab,f);
  P = epg_X_grad(P,1);
  P = epg_X_rf(P,abs(flipangle(ech)),angle(flipangle(ech)));   
  P = epg_X_grad(P,1);
  P = epg_X_relax(P,esp/2,T1,T2,ka,deltab,f);
  s(1:2,ech) = P([1 4],1);  	% Signal is F0 state.
  s(3,ech) =  sum(P([1 4],1));
end

end