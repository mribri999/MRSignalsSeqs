
% Lecture 8, Example 02
% 
% Calculation of basic signal in a spin-echo sequence using EPG

T1=2;	% sec
T2=0.1	% sec
TR=1;	% sec
TE=0.05;% sec
N = 2	% cycles (Nth cycle signal saved)

Q = epg_m0();	% Start at equilibrium

for n=1:N	% Repeat twice, to get to steady-state

  Q = epg_rf(Q,pi/2,pi/2);	% Flip 90 about y
  Q = epg_grelax(Q,T1,T2,TE/2,0,0,1);	% Gradient to TE/2
  Q = epg_rf(Q,pi,0);	% Flip 180 about x
  Q = epg_grelax(Q,T1,T2,TE/2,0,0,1);	% Gradient to TE

  Qecho = Q;	% Save spin-echo signal

  Q = epg_grelax(Q,T1,T2,TR-TE,0,0,1);	% TE to TR.  
					% Set last 1 to 0 to not dephase

end;

tt = sprintf('Spin-Echo signal on TR %d is %g',N,Qecho(1)); disp(tt);






