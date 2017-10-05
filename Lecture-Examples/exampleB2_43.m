% Stimulated Echo Sequence (Bloch vs EPG) example
%
%	Slides contain matrix examples
%
figure(1);

FZ = [0;0;1]			% Equilibrium
epg_show(FZ,0,1,23); disp('0 - Press Enter'); drawnow; pause;
FZ = epg_rf(FZ,pi/2,pi/2)	% After 90y
epg_show(FZ,0,1,23);  disp('1 - Press Enter'); drawnow; pause;
FZ = epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
epg_show(FZ,0,1,23); disp('2 - Press Enter'); drawnow; pause;
M = epg_FZ2spins(FZ,12);	% Check!
FZ = epg_rf(FZ,pi/2,pi/2)	% After 90y
epg_show(FZ,0,1,23); disp('3 - Press Enter'); drawnow; pause;
FZ = epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
epg_show(FZ,0,1,23); disp('4 - Press Enter'); drawnow; pause;
FZ = epg_rf(FZ,pi/2,pi/2)	% After 90y
epg_show(FZ,0,1,23); disp('5 - Press Enter'); drawnow; pause;
FZ = epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
epg_show(FZ,0,1,23); disp('6 - Press Enter'); drawnow; pause;











