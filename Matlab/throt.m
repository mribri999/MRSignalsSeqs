%	function [M] = throt(alpha,phi)
%
%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the axis
%	defined by y=x*tan(phi), by alpha degrees.
%	The rotation is left-handed.
%
%	phi is in degrees too.
%
%	See also: tprot
%

% ======================== CVS Log Messages ========================
% $Log: throt.m,v $
% Revision 1.2  2002/03/28 00:50:11  bah
% Added log to source file
%
%
%
% ================================================================== 


function [M] = throt(alpha,phi)

ca = cos(pi/180*alpha);	% cosine of tip alpha
sa = sin(pi/180*alpha);	% sine of tip
cp = cos(pi/180*phi  ); % cosine of phi
sp = sin(pi/180*phi  ); % sine of phi


M = [cp*cp+sp*sp*ca cp*sp*(1-ca) -sp*sa;
     cp*sp-sp*cp*ca sp*sp+cp*cp*ca cp*sa;
     sa*sp -sa*cp ca];








