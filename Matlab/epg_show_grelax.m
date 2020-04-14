%function [FpFmZ] = epg_show_grelax(FpFmZ,noadd,T1,T2,T,Nanim,showtwists,shiftfs)
%
%	Propagate EPG states through a "unit" gradient, with relaxation
%	and animate.
%	NOTE:  This is a bit contrived, as the magnetization is not
%	really all in F0 (say) as it dephases.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%               T1,T2 = Relaxation times (same as T)
%               T = Time interval (same as T1,T2)
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
function [FpFmZ] = epg_show_grelax(FpFmZ,noadd,T1,T2,T,Nanim,showtwists,shiftfs)

if (nargin < 1 || length(FpFmZ)<1) FpFmZ = [.3 .5 ;.3 .5 ; .2 .2 ]; end;
if (nargin < 2 || length(noadd)<1) noadd=0; end;	% Add by default.  
if (nargin < 3 || length(T1)<1) T1 = 1; end;
if (nargin < 4 || length(T2)<1) T2 = 0.2; end;
if (nargin < 5 || length(T)<1) T = .2; end;
if (nargin < 6 || length(Nanim)<1) Nanim=48; end;
if (nargin < 7 || length(showtwists)<1) showtwists=0; end;
if (nargin < 8 || length(shiftfs)<1) shiftfs=0; end;

[m,n] = size(FpFmZ);
if (m<3) FpFmZ(3,1)=0; end;

% -- Added for animation... 
FZ = FpFmZ;
for k=1:Nanim
  FZ = epg_relax(FZ,T1,T2,T/Nanim); 
  epg_show(FZ(1:m,:),k/Nanim,[],[],showtwists,shiftfs);
  drawnow;
end;


% -- Actually DO epg_grad and relaxation!
FpFmZ = epg_grad(FpFmZ,noadd);
FpFmZ = epg_relax(FpFmZ,T1,T2,T);

