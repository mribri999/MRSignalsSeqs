%
%	Just play with animation functions a bit more, showing 
%	all states.

FZ = [0;0;1];	% M0 (along Z)
FZ = epg_show_rf(FZ,-pi/2,0);		% Along x!
disp('Resize plot if you like, and press Enter!');
pause;
FZ = epg_show_grad(FZ);		% Gradient
FZ = epg_trim(FZ,0.05);
FZ = epg_show_grad(FZ);		% Another Gradient
FZ = epg_trim(FZ,0.05);

FZ = epg_show_rf(FZ,-pi,pi/2);		% Along y!

FZ = epg_show_grad(FZ);		% Gradient
FZ = epg_trim(FZ,0.05);
FZ = epg_show_grad(FZ);		% Another Gradient
FZ = epg_trim(FZ,0.05);

 
