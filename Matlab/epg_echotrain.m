%	
%	epg_echotrain.m - Example script
%	Simulate spin echo trains with different refocusing pulses.
%
% 

s1 = epg_cpmg(pi,30,1000,200,10);	% IDEAL, CPMG
s2 = epg_cpmg(pi*i,30,1000,200,10);	% IDEAL, non CPMG
s3 = epg_cpmg(2/3*pi,30,1000,200,10);	% 120 CPMG
flips = 2/3*pi*ones(1,30);		% Constant flip angles
s4 = epg_cpmg(flips,30,1000,200,10);	% 120 without 1st pulse correction.
s5 = epg_cpmg(2/3*pi*i,30,1000,200,10);	% 120 non CPMG


t = [1:30]*10;	% Echo times.

plot(t,abs([s1(:) s2(:) s3(:) s4(:) s5(:)]));
legend('180 CPMG','180 Non-CPMG','120 CPMG','120 CPMG no corr','120 Non-CPMG');





