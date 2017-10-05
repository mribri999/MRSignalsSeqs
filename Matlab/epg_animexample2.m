%
%	Example of EPG animations (Gradient-spoiled sequence)
%
%	Some parameters:  (Remove % as needed!)

flipa= 60*i*pi/180;
T1 = 500;
T2 = 100;
TR = 5;
% TE = 0;

Nframes = 5;
scale = 0.5;
Nspins = 30;

FZ = [0;0;1];	% Equilibrium

% Let Steady State form

for k=1:200
  FZ = epg_rf(FZ,abs(flipa),angle(flipa));	% alpha pulse.
  FZ = epg_grelax(FZ,T1,T2,TR,0,0,1);
end;



for k=1:Nframes
	% 90 degree pulse.
  FZ = epg_animrf(FZ,abs(flipa),angle(flipa),Nspins,scale);

  title('Echo Signal');
  pause(1);
  FZ = epg_animgrad(FZ,Nspins,scale);			% Spoiler
end;




