
[acq, sys, Gx, Gy, Gz, RF] = Rad229_Slice_Select_Demo;

zlocs = [-20:20];	 % mm
dT = .01;		 % ms

rf = 1000000*RF.B1;	 % Convert T to uT
G = Gz.G * 1000;	 % Convert T/m to mT/m

figure;

%for n=length(rf) 	% Just show end point.
for n=1:3:length(rf)	% Repeat for all points.
  [mx,my,mz] = blochsim(rf(1:n),G(1:n),dT,0,zlocs);
  plotc(zlocs,mx+i*my);
  lplot('Z location','Mxy','Mxy after Slice Selection');
  a = axis; axis([a(1:2) -1 1]);	% Show magnitude -1 to 1
  pause(0.1);
end;

% Plot portions along the way...?

