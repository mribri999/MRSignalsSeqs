% Lecture 9, Example 03
%
% Signals vs flip angled for bSSFP, gradient-spoiled, RF-spoiled.

TR = 5;         % ms
TE = 0; 	% ms
flips = [[0:0.1:5] [6:60]];      % degrees
T1 = 1000;      % ms
T2 = 200;       % ms

sbssfp0 = zeros(length(flips));
sbssfp180 = sbssfp0;
sgre= sbssfp0;
srfsp= sbssfp0;

for f=1:length(flips)
  sbssfp0(f) = bssfp(flips(f),TR,TE,T1,T2,0);
  sbssfp180(f) = bssfp(flips(f),TR,TE,T1,T2,500/TR);
  sgre(f) = gresignal(T1,T2,TE,TR,flips(f));
  srfsp(f) = exrecsignal(T1,T2,TE,TR,flips(f));
end;

clf;
set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width
figure(1);
plot(flips,abs(sbssfp0),flips,abs(sbssfp180),flips,abs(sgre),flips,abs(srfsp));
lplot('Flip Angle (deg)','Signal Magnitude','Gradient-Echo Signals vs Flip Angle');
legend('bSSFP (0^\circ)','bSSFP (180^\circ)','Gradient-Spoiled','RF-spoiled');


