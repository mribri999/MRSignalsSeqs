%
%function [FpFmZ] = epg_mgrad(FpFmZ,noadd)
%	Propagate EPG states through a *NEGATIVE* "unit" gradient.
%	(States move the opposite direction to epg_grad!)
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%
%	OUTPUT:
%		Updated FpFmZ state.
%
%       SEE ALSO:
%               epg_grad, epg_grelax
%
%       B.Hargreaves.
function [FpFmZ] = epg_mgrad(FpFmZ,noadd)

if (nargin < 2) noadd=0; end;	% Add by default.  

% Gradient does not affect the Z states.

if (noadd==0)
  FpFmZ = [FpFmZ [0;0;0]];	% Add higher dephased state.
end;

FpFmZ(2,:) = circshift(FpFmZ(2,:),[0 1]);	% Shift Fm states.
FpFmZ(1,:) = circshift(FpFmZ(1,:),[0 -1]);	% Shift Fp states.
FpFmZ(1,end)=0;					% Zero highest Fp state.
FpFmZ(2,1) = conj(FpFmZ(1,1));			% Fill in lowest Fm state.



