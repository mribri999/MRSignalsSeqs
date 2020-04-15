%
%	Just play with animation functions a bit more, showing 
%	all states.  This is a sequence with 2 spin-echoes
pauseN = 12;	% pause at each point.

FZ = epg_m0(2);	% M0 (along Z)
FZ = epg_show_rf(FZ,pi/2,pi/2);		% M Along x!
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_rf(FZ,pi,0);		% Along y!
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_rf(FZ,pi,0);		% Along y!
for k=1:pauseN; epg_saveframe; end;

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
for k=1:pauseN; epg_saveframe; end;

