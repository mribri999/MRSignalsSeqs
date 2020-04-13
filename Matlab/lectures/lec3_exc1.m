%	Lecture 3, Exercise 1
%
%	Aim is to plot the flip angle vs T2 for a rectangular
%	RF pulse of duration 0.5ms, and amplitude ~12uT.  Nominal
%	flip angle is 90 degrees.
%
%

Nrf = 1000;			% points
Trf = 0.0005;			% sec
gbar = 42.58;			% kHz/mT
Arf = 0.25/gbar/(1000*Trf);	% mT

T1 = 1.1;			% sec
T2s = logspace(-5,0,40);	% T2 values to test

peakflip = Trf * gbar * Arf;	% Rotation in rotations
Mxy = zeros(length(T2s));	% Allocate

for t=1:length(T2s)

end;

figure(1);
semilogx(T2s,Mxy,T2s,alph/(pi/2));
lplot('T2 (sec)','M_{xy}/M_0','Excitation vs T_2');
legend('M_{xy}/M_0','Flip Angle / 90^\circ');

	




