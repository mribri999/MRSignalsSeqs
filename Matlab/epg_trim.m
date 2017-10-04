function [FpFmZ] = epg_trim(FpFmZ,thres)
%
%	Trim higher-order states to N if the sum of absolute
%	values of F+, F- and Z of each order is less than
%	thres for all orders n>N.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		thres = threshold for trimming.
%       OUTPUT:
%               Updated FpFmZ state.
%
%       B.Hargreaves.


f = find(sum(abs(FpFmZ))>=thres);
fn = max(f);
FpFmZ = FpFmZ(:,1:fn);

