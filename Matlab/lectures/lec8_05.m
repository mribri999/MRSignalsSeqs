
% Lecture 8, Example 05
% 
% Try to stabilize the 60-degree refocusing flip angle signal

clf;
flips = pi/180*[120 60*ones(1,23)];
set(0,'defaultAxesFontSize',14);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width

s = epg_cpmg(flips);



