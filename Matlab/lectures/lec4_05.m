
% Lecture 4, Example 05
% 
% A stimulated echo - simple simulation

Q = epg_m0(5);
figure(1);
Q = epg_rf(Q,pi/2,0); figure(1); epg_show(Q);
Q = epg_grad(Q,1); figure(2); epg_show(Q);
Q = epg_rf(Q,pi/2,0); figure(3); epg_show(Q);
Q = epg_grad(Q,1); figure(4); epg_show(Q);
Q = epg_grad(Q,1); figure(5); epg_show(Q);
Q = epg_grad(Q,1); figure(6); epg_show(Q);
Q = epg_rf(Q,pi/2,0); figure(7); epg_show(Q);
Q = epg_grad(Q,1); figure(8); epg_show(Q);


