%	function [M] = xrot(angle)
%
%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the x axis
%	by angle degrees, using a left-handed rotation.
%

% ======================== CVS Log Messages ========================
% $Log: xrot.m,v $
% Revision 1.2  2002/03/28 00:50:12  bah
% Added log to source file
%
%
%
% ================================================================== 

function [M] = xrot(angle)

c = cos(pi*angle/180);
s = sin(pi*angle/180);

M = [1 0 0; 0 c s; 0 -s c];



