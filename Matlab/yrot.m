%	function [M] = yrot(angle)
%
%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the y axis
%	by angle degrees, using a left-handed rotation.
%

% ======================== CVS Log Messages ========================
% $Log: yrot.m,v $
% Revision 1.2  2002/03/28 00:50:12  bah
% Added log to source file
%
%
%
% ================================================================== 

function [M] = yrot(angle)

c = cos(pi*angle/180);
s = sin(pi*angle/180);

M = [c 0 -s; 0 1 0; s 0 c];




