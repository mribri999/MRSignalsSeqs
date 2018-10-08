%
%function [FpFmZ] = epg_grad_show(FpFmZ,noadd)
%	Propagate EPG states through a "unit" gradient, and animate as it
%	goes.  NOTE:  This is a bit contrived, as the magnetization is not
%	really all in F0 (say) as it dephases.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%
%	OUTPUT:
%		Updated FpFmZ state.
%
%       SEE ALSO:
%               epg_grad, epg_grelax
%
%       B.Hargreaves.
function [FpFmZ] = epg_grad(FpFmZ,noadd)

% -- Added for animation... 
Nanim=16;	% Animation frames
for k=1:Nanim
  epg_show(FpFmZ,k/Nanim);
  drawnow;
end;


% -- Actually DO epg_grad!
if (nargin < 2) noadd=0; end;	% Add by default.  
FpFmZ = epg_grad(FpFmZ,noadd);
epg_show(FpFmZ);
drawnow;

