%	function [M] = zrot(angle)
%
%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the z axis
%	by angle degrees, using a left-handed rotaton.
%

% ======================== CVS Log Messages ========================
% $Log: zrot.m,v $
% Revision 1.2  2002/03/28 00:50:12  bah
% Added log to source file
%
%
%
% ================================================================== 


function [M] = zrot(angle)

c = cos(pi*angle/180);
s = sin(pi*angle/180);

M = [c s 0; -s c 0; 0 0 1];


