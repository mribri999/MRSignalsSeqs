% Lecture 4, Example 8

% Calculate steady state from RF + gradient

N = 6;
Q = epg_m0(N);		% Equilibrium magnetization, use N states.

for k=1:1000;
  Q = epg_grelax(Q,1,.3,.005,0,0,1,1);	% Gradient and relaxation
  %epg_show(Q(:,1:6));
  %drawnow(); pause;
  Q = epg_rf(Q,pi/6,pi/2);		% Apply RF rotation
end;

epg_show(Q(:,1:4));		% Show result
% epg_showstate(Q);		% Show 3D magnetization

