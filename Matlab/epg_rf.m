%
%function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)
%	Propagate EPG states through an RF rotation of 
%	alpha, with phase phi (both radians).
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%
%       OUTPUT:
%               FpFmZ = Updated FpFmZ state.
%		RR = RF rotation matrix (3x3).
%

%	SEE ALSO:

%		epg_grad, epg_grelax
%
%	B.Hargreaves.
%
function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)

if (nargin < 1) FpFmZ = [0;0;1];	end;	% Default is M=M0
if (nargin < 2) alpha = pi/2;		end;	% Default is pi/2
if (nargin < 3) phi = 0;		
  tt = sprintf('epg_rf defaulting to M0 with alpha=%g, phi=%g',alpha,phi);
  disp(tt);
end;	% Default is pi/2

% -- From Weigel at al, JMR 205(2010)276-285, Eq. 8.

if (nargin < 3) phi = pi/2; end;
if length(FpFmZ)<3 FpFmZ = [0;0;1]; end;

if (abs(alpha)>2*pi) warning('epg_rf:  Flip angle should be in radians!'); end;


RR = [(cos(alpha/2))^2 exp(2*i*phi)*(sin(alpha/2))^2 -i*exp(i*phi)*sin(alpha);
      exp(-2*i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 i*exp(-i*phi)*sin(alpha);
      -i/2*exp(-i*phi)*sin(alpha) i/2*exp(i*phi)*sin(alpha)      cos(alpha)];


FpFmZ = RR * FpFmZ;


