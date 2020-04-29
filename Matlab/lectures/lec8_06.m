
% Lecture 8, Example 06
% 
% Hyperecho example.

% -- These lines setup a short hyperecho sequence:
flips = [1 1.3 1.2 1.3] .* exp(pi*i*[0.62 -.28 .28 .25]);	% Random
%flips = [2 1 1 1 1 1];		% Just reduced (~60 degree) train



flips=[flips pi -fliplr(conj(flips))];         % Hyperecho refocusing

clf;
set(0,'defaultAxesFontSize',14);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width
s = epg_cpmg(flips);	



