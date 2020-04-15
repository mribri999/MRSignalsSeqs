% Lecture 4, Example 01
% 
% Simple image reconstruction in matlab
%
[Q,R] = epg_rf(epg_m0(), pi,0); disp('Rotation about Mx:'); disp(R);
[Q,R] = epg_rf(epg_m0(), pi,pi/2); disp('Rotation about My:'); disp(R);


