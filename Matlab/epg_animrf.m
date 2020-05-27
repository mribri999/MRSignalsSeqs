
function FZ = epg_animrf(FZ,flip,phase,N,scale,showorder)

%
%	Animate an EPG state through an RF pulse
%

if (nargin < 4) N = 24; end;
if (nargin < 5) scale = 1; end;
if (nargin < 6) showorder=-1; end;

Nf = 20;

for p=0:Nf

  if (showorder==-1)
    epg_showstate(FZ,0,scale,N);
  else
    epg_showorder(FZ(:,showorder+1),showorder,0);
  end;

  if (p<Nf)
    FZ = epg_rf(FZ,flip/Nf,phase);
  end;

  drawnow;
  pause(0.05);
  
  if (exist('epg_saveframe'))
    epg_saveframe;  % Save frame, if globals 'framenum' and 'filestem' exist
  end;


end;



