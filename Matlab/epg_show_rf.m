%function [FpFmZ,RR] = epg_show_rf(FpFmZ,alpha,phi,Nanim,showtwists)
%
%	Propagate EPG states through an RF rotation of 
%	alpha, with phase phi (both radians), plotting frames.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		alpha, phi = flip angle and phase.
%               Nanim = Number of frames
%               showtwists = 1 to show twist, 0 to have all spins from origin
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
function [FpFmZ,RR] = epg_show_rf(FpFmZ,alpha,phi,Nanim,showtwists)

if (nargin < 1 || length(FpFmZ)<1) FpFmZ = [0;0;1]; end;
if (nargin < 3 || length(phi)<1) phi = pi/2; end;
if (nargin < 2 || length(alpha)<1) alpha = pi/2; end;		
if (nargin < 4 || length(Nanim)<1) Nanim=32; end;
if (nargin < 5 || length(showtwists)<1) showtwists=0; end;

[m,n] = size(FpFmZ);
if (m<3) FpFmZ(3,1)=0; end;

% -- From Weigel at al, JMR 205(2010)276-285, Eq. 8.

if (abs(alpha)>2*pi) warning('epg_rf:  Flip angle should be in radians!'); end;

RR = [(cos(alpha/2))^2 exp(2*i*phi)*(sin(alpha/2))^2 -i*exp(i*phi)*sin(alpha);
      exp(-2*i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 i*exp(-i*phi)*sin(alpha);
      -i/2*exp(-i*phi)*sin(alpha) i/2*exp(i*phi)*sin(alpha)      cos(alpha)];

alpha = alpha/Nanim;
RRa = [(cos(alpha/2))^2 exp(2*i*phi)*(sin(alpha/2))^2 -i*exp(i*phi)*sin(alpha);
      exp(-2*i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 i*exp(-i*phi)*sin(alpha);
      -i/2*exp(-i*phi)*sin(alpha) i/2*exp(i*phi)*sin(alpha)      cos(alpha)];

FZ = FpFmZ;

% -- Animate.
for k=1:Nanim
  FZ = RRa * FZ;
  epg_show(FZ(1:m,:),[],[],[],showtwists);
  drawnow;
end;

FpFmZ = RR * FpFmZ;


