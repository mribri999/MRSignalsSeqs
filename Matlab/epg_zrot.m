%function [FpFmZ] = epg_zrot(FpFmZ,rotangle)
%
%	Rotate spin state about Mz axis by rotangle radians.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		rotangle = angle by which to rotate (radians)
%
%	OUTPUT:
%		Updated FpFmZ state.
%
%       SEE ALSO:
%               epg_grad, epg_grelax
%
%       B.Hargreaves.
%
function [FpFmZ] = epg_zrot(FpFmZ,rotangle)

if (nargin < 2) rotangle = 0; end;
phasor = exp(i*rotangle);

FpFmZ = diag([phasor conj(phasor) 1]) * FpFmZ;

