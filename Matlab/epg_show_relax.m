%function [FpFmZ] = epg_show_relax(FpFmZ,T1,T2,T,Nanim,showtwists)
%
%	Propagate EPG states through relaxation given by T1,T2,T
%	plotting frames.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%               T1,T2 = Relaxation times (same as T)
%               T = Time interval (same as T1,T2)
%               Nanim = Number of frames
%               showtwists = 1 to show twist, 0 to have all spins from origin
%
%       OUTPUT:
%               FpFmZ = Updated FpFmZ state.
%
%	SEE ALSO:
%		epg_grad, epg_grelax
%
%	B.Hargreaves.
%
function [FpFmZ] = epg_show_relax(FpFmZ,T1,T2,T,Nanim,showtwists)

if (nargin < 1 || length(FpFmZ)<1) FpFmZ = [.3 .5 .25;.3 .5 .25; .2 .2 .2]; end;
if (nargin < 2 || length(T1)<1) T1 = 1; end;		
if (nargin < 3 || length(T2)<1) T2 = 0.2; end;
if (nargin < 4 || length(T)<1) T = 1; end;
if (nargin < 5 || length(Nanim)<1) Nanim=16; end;
if (nargin < 6 || length(showtwists)<1) showtwists=0; end;

T = T/Nanim;

% -- Animate.
for k=1:Nanim
  FpFmZ = epg_relax(FpFmZ,T1,T2,T);
  epg_show(FpFmZ,[],[],[],showtwists);
  drawnow;
end;


