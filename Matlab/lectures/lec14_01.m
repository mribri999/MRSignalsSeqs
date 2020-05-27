% Lecture 14, Example 01
%
% Animate saturation pulse with recovery.

global filestem;
global framenum;
filestem='/Users/brian/tmp/sat';
framenum=0;


FZ = epg_m0;

% == Animate RF
Nrf = 24;
for k=0:Nrf 
  epg_showorder(FZ,0);
  if (k<Nrf) FZ = epg_rf(FZ,pi/2/Nrf,pi/2); end;
  subplot(2,2,4); lplot('','',''); title(''); axis off;
  drawnow;
  epg_saveframe;
end;

% == Animate Gradient
Ng = 24;
for k=1:Ng
  epg_showorder(FZ,0,k/Ng,23);
  subplot(2,2,4); lplot('','',''); title(''); axis off;
  drawnow;
  epg_saveframe;
end;

% == Animate Relaxation
Nr = 100;
for k=1:Nr
  FZ = epg_grelax(FZ,300,80,6,0,0,0,1);		% Relaxation
  epg_showorder(FZ,0,1,23);
  subplot(2,2,4); lplot('','',''); title(''); axis off;
  drawnow;
  epg_saveframe;
end;



  
