%
%	Just play with animation functions a bit more, showing 
%	all states.  This is a sequence with 2 spin-echoes
pauseN = 12;	% pause at each point.

refocflip = 130/180*pi;

FZ = epg_m0(3);	% M0 (along Z)
FZ = epg_show_rf(FZ,pi/2,pi/2);		% M Along x!

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient

FZ = epg_show_rf(FZ,refocflip,0);		% Along y!

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient

FZ = epg_show_rf(FZ,refocflip,0);		% Along y!

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient

FZ = epg_show_rf(FZ,refocflip,0);		% Along y!

FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient
FZ = epg_show_grelax(FZ,1,195,95,10,36,1,1);		% Gradient

FZ = epg_show_rf(FZ,refocflip,0);		% Along y!

