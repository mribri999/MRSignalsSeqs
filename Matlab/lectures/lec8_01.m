
% Lecture 8, Example 01
% 
% Simulate a simple 90-gradient-RF-gradient to a spin echo, using EPG

refoc_angle = pi;

figure(1); clf;

subplot(3,2,1);
Q = epg_m0();   % Start at equilibrium
epg_showstate(Q); title('Equilibrium'); 

subplot(3,2,2);
Q = epg_rf(Q,pi/2,pi/2);      % Flip 90 about y
epg_showstate(Q); title('After Excitation'); 

subplot(3,2,3);
Q = epg_grad(Q);
epg_showstate(Q); title('After Gradient'); 

subplot(3,2,4);
Q = epg_rf(Q,refoc_angle,0);   % Flip 180 about x
epg_showstate(Q); title('After Refocusing'); 

subplot(3,2,5);
Q = epg_grad(Q);
epg_showstate(Q); title('Spin Echo'); 






