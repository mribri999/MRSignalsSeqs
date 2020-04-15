% Lecture 4, Example 02
% 
% Simple image reconstruction in matlab
%
disp('-- Spin Echo with 180x pulses --')
Q = epg_m0(2); disp('Equilibrium: '); disp(Q)
T1 = 195; T2=95; T=10;     % Chosen for E1=0.95, E2=0.9
Q = epg_relax(Q,T1,T2,T);
Q = epg_rf(Q,pi/2,pi/2); disp('Excitation: '); disp(Q)
Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('Pre-180: '); disp(Q)
Q = epg_rf(Q,pi,0); disp('Post-180'); disp(Q)
Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('1st Spin Echo: '); disp(Q)
Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('Pre-180: '); disp(Q)
Q = epg_rf(Q,pi,0); disp('Post-180'); disp(Q)
Q = epg_relax(epg_grad(Q,1),T1,T2,T); disp('2nd Spin Echo: '); disp(Q)

