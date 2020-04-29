
% Lecture 8, Example 07
% 
% Two spin-echoes, as in the lecture example, with/without CPMG and
% with/without stabilization.

refoc_angle = 120/180*pi;
stabilize=1;
cpmg=0;
pausit=1;

figure(1); clf;

Q0 = epg_m0();   % Start at equilibrium
epg_show(Q0); title('Equilibrium'); disp('press enter'); pause;
if (pausit) disp('press enter'); pause; end;

if (cpmg)
  Q0= epg_rf(Q0,pi/2,pi/2)      % Flip 90 about y
else
  Q0= epg_rf(Q0,pi/2,0)      % Flip 90 about x 
end;
epg_show(Q0); title('After Excitation'); 
if (pausit) disp('press enter'); pause; end;

Q1 = epg_grad(Q0)
epg_show(Q1); title('After Gradient'); 
if (pausit) disp('press enter'); pause; end;

if (stabilize==0)
  Q2 = epg_rf(Q1,refoc_angle,0)   % Refocus about x
else
  Q2 = epg_rf(Q1,(pi+refoc_angle)/2,0)   % Refocus about x
end;

epg_show(Q2); title('After Refocusing'); 
if (pausit) disp('press enter'); pause; end;

Q3 = epg_grad(Q2)
epg_show(Q3); title('Spin Echo'); 
if (pausit) disp('press enter'); pause; end;
%Q3(1)=0;	% Zero F0 state
Q3(3,2)=0;	% Zero Z1 state.

Q4 = epg_grad(Q3)
epg_show(Q4); title('before next refocusing'); 
if (pausit) disp('press enter'); pause; end;

Q5 = epg_rf(Q4,refoc_angle,0)   % Refocus about x
epg_show(Q5); title('After Refocusing'); 
if (pausit) disp('press enter'); pause; end;

Q6 = epg_grad(Q5)
epg_show(Q6); title('before next refocusing'); 
if (pausit) disp('press enter'); pause; end;




