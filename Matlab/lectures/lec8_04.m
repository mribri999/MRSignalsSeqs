
% Lecture 8, Example 04
% 
% First demonstration of epg_cpmg.m for constant refocusing pulses with stabilization and CPMG

flips = 20*[1:9];	% Repeat for different refocusing angles
etl = 24;
set(0,'defaultAxesFontSize',14);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width
fliplist = {};

figure(1); clf
for k=1:length(flips)
  s = epg_cpmg(flips(k)*pi/180,etl);	% Will plot each time, just keep last.
  sigs(:,k) = s(:);
  fliplist = {fliplist{:}, sprintf('%d',flips(k))};
end;

figure(2);
plot([1:etl],abs(sigs));
lplot('Echo#','Signal','CPMG Signal vs Refocusing Angle');
legend(fliplist);


