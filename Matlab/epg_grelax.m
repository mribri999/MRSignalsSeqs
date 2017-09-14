%function [FpFmZ,EE,BV] = epg_grelax(FpFmZ,T1,T2,T,kg,D,Gon,noadd)
%
%	Propagate EPG states through a period of relaxation, and
%	diffusion over an interval T, with or without a gradient.
%	Leave last 3 blank to exclude diffusion effects.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		T1,T2 = Relaxation times (s)
%		T = Time interval (s)
%	    (Optional inputs follow)
%		kg = k-space traversal due to gradient (rad/m) for diffusion
%		D = Diffusion coefficient (m^2/s)
%		Gon = 0 if no gradient on, 1 if gradient on
%			(gradient will advance states at the end.)
%		noadd=1 to not add higher-order states - see epg_grad.m
%
%	OUTPUT:
%		FpFmZ = updated F+, F- and Z states.
%		EE = decay matrix, 3x3 = diag([E2 E2 E1]);
%		BV = b-value matrix, 3xN (see FpFmZ) of attenuations.
%
%       SEE ALSO:
%               epg_grad, epg_rf
%
%       B.Hargreaves.

function [FpFmZ,EE,BV] = epg_grelax(FpFmZ,T1,T2,T,kg,D,Gon,noadd)

if (nargin < 5) kg = 0; D = 0; end;	% Default ignore diffusion
if (nargin < 8) noadd=0; end;	% Default is to add states.
if (nargin < 7) Gon = 1; end;	% Default is Gon.

if (nargin > 1)			% Skip relaxation if only one argument
  E2 = exp(-T/T2);
  E1 = exp(-T/T1);

  EE = diag([E2 E2 E1]);	% Decay of states due to relaxation alone.
  RR = [1-E1];			% Mz Recovery, affects only Z0 state, as 
				% recovered magnetization is not dephased.


  FpFmZ = EE * FpFmZ;		% Apply Relaxation
  FpFmZ(3,1) = FpFmZ(3,1)+RR;	% Recovery  ( here applied before diffusion,
				% but could be after or split.)
end;


if (nargin > 4)			% Model Diffusion Effects

  Findex = 0:length(FpFmZ(1,:))-1;	% index of states, 0...N-1
  bvalZ = ((Findex)*kg).^2*T;		% diffusion  for Z states, assumes that
					% the Z-state has to be refocused, so
					% this models "time between gradients"
	
	% For F states, the following models the additional diffusion time
	% (Findex) and the fact that the state will change if the gradient is
	% on (0.5*Gon), then the additional diffusion *during* the gradient, 
	% ... Gon*kg^2/12 term.

  bvalp = ((( Findex+.5*Gon)*kg).^2+Gon*kg^2/12)*T;	% for F+ states
  bvalm = (((-Findex+.5*Gon)*kg).^2+Gon*kg^2/12)*T;	% for F- states

		

  FpFmZ(1,:) = FpFmZ(1,:) .* exp(-bvalp*D);	% diffusion on F+ states
  FpFmZ(2,:) = FpFmZ(2,:) .* exp(-bvalm*D);	% diffusion on F- states
  FpFmZ(3,:) = FpFmZ(3,:) .* exp(-bvalZ*D);	% diffusion of Z states.
 
  BV = [bvalp; bvalm; bvalZ];	% For output. 
end;

if (Gon==1)
  if (kg >= 0)
    FpFmZ = epg_grad(FpFmZ,noadd);	% Advance states.
  else 
    FpFmZ = epg_mgrad(FpFmZ,noadd);	% Advance states by negative gradient.
  end;
end;


