
function FZ = epg_animgrad(FZ,N,scale)
%function FZ = epg_animgrad(FZ,N,scale)
%
%	Animate an EPG state to the next F state.
%
%	N = #spins
%	scale = axis scale multiplier
%

if (nargin < 2) N = 24; end;
if (nargin < 3) scale = 1; end;

Nf = 20;

for p=1:Nf

  epg_showstate(FZ,p/Nf,scale,N);
  drawnow;
  pause(0.05);

  if (exist('epg_saveframe'))
    epg_saveframe;  % Save frame, if globals 'framenum' and 'filestem' exist 
  end;

end;

FZ = epg_grad(FZ);


