%function [FpFmZ] = epg_show_grad(FpFmZ,noadd,Nanim,showtwists)
%
%	Propagate EPG states through a "unit" gradient, and animate as it
%	goes.  NOTE:  This is a bit contrived, as the magnetization is not
%	really all in F0 (say) as it dephases.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%		Nanim = Number of frames
%		showtwists = 1 to show twist, 0 to have all spins from origin
%
%	OUTPUT:
%		Updated FpFmZ state.
%
%       SEE ALSO:
%               epg_grad
%
%       B.Hargreaves.
%
function [FpFmZ] = epg_show_grad(FpFmZ,noadd,Nanim,showtwists)

if (nargin < 2) noadd=0; end;	% Add by default.  
if (nargin < 3 || length(Nanim)<1) Nanim=32; end;
if (nargin < 4 || length(showtwists)<1) showtwists=0; end;

% -- Added for animation... 
Nanim=16;	% Animation frames
for k=1:Nanim
  epg_show(FpFmZ,k/Nanim,[],[],showtwists);
  drawnow;
end;


% -- Actually DO epg_grad!
FpFmZ = epg_grad(FpFmZ,noadd);
epg_show(FpFmZ);
drawnow;

