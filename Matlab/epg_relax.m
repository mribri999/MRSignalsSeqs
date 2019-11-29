%function [FpFmZ,EE] = epg_relax(FpFmZ,T1,T2,T)
% 
%	Propagate EPG states through a period of relaxation over
%	an interval T.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		T1,T2 = Relaxation times (same as T)
%		T = Time interval (same as T1,T2)
%
%	OUTPUT:
%		FpFmZ = updated F+, F- and Z states.
%		EE = decay matrix, 3x3 = diag([E2 E2 E1]);
%
%       SEE ALSO:
%               epg_grad, epg_rf
%
%       B.Hargreaves.

function [FpFmZ,EE] = epg_gt(FpFmZ,T1,T2,T)

if(T1 < 0 || T2 < 0 || T < 0)
    warning('negative values for time')
end

E2 = exp(-T/T2);
E1 = exp(-T/T1);

EE = diag([E2 E2 E1]);	% Decay of states due to relaxation alone.
RR = [1-E1];			% Mz Recovery, affects only Z0 state, as 
				% recovered magnetization is not dephased.


FpFmZ = EE * FpFmZ;		% Apply Relaxation
FpFmZ(3,1) = FpFmZ(3,1)+RR;	% Recovery  


