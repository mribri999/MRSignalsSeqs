
% Lecture 8, Example 03
% 
% Simulate a simple 90-gradient-RF-gradient to a spin echo, using EPG
% Same as example 1, but showing all EPG states.

refoc_angle = 120/180*pi;

figure(1); clf;

Q = epg_m0();   % Start at equilibrium
epg_show(Q); title('Equilibrium'); disp('press enter'); pause;

Q = epg_rf(Q,pi/2,pi/2);      % Flip 90 about y
epg_show(Q); title('After Excitation'); disp('press enter'); pause;

Q = epg_grad(Q);
epg_show(Q); title('After Gradient'); disp('press enter'); pause;

Q = epg_rf(Q,refoc_angle,0);   % Refocus about x
epg_show(Q); title('After Refocusing'); disp('press enter'); pause;

Q = epg_grad(Q);
epg_show(Q); title('Spin Echo'); disp('press enter'); pause;






