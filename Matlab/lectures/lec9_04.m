% Lecture 9, Example 04
%
% EPG simulation of Signals vs flip angled for bSSFP, gradient-spoiled, RF-spoiled.
% Notice that this calls ONE function that simulates any sequence iteratively
% to the steady state.  It should clarify the sequence steps for each 
% gradient-echo variant.

TR = 5;         % ms
TE = 0;         % ms
flips = [[0:0.1:5] [6:60]]*pi/180;      % radians
T1 = 1000;      % ms
T2 = 200;       % ms

sbssfp0 = zeros(length(flips));
sbssfp180 = sbssfp0;
sgre= sbssfp0;
srfsp= sbssfp0;

for f=1:length(flips)
  sig = epg_gradecho(flips(f),T1,T2,TR,TE,0);
  sbssfp0(f) = sig(end);
  sig = epg_gradecho(flips(f),T1,T2,TR,TE);
  sbssfp180(f) = sig(end);
  sig = epg_gradecho(flips(f),T1,T2,TR,TE,0,0,1);
  sgre(f) = sig(end); 
  sig = epg_gradecho(flips(f),T1,T2,TR,TE,0,117/180*pi,1);
  srfsp(f) = sig(end);
end;

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width
figure;
plot(flips,abs(sbssfp0),flips,abs(sbssfp180),flips,abs(sgre),flips,abs(srfsp));
lplot('Flip Angle (deg)','Signal Magnitude','EPG-Simulated Gradient-Echo Signals vs Flip Angle');
legend('bSSFP (0^\circ)','bSSFP (180^\circ)','Gradient-Spoiled','RF-spoiled');


