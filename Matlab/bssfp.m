% function s = bssfp(a,TR,TE,T1,T2,df)
%
%	Calculate the bSSFP signal as a function of parameters
%
%	a = flip angle in degrees
%	Times in ms, df in Hz

function s = bssfp(flip,TR,TE,T1,T2,df)

E1 = exp(-TR/T1);
E2 = exp(-TR/T2);

flip = flip*pi/180;
phi = 2*pi*df*TR/1000;	% Precession over TR (rad)

sf = sin(flip);
cf = cos(flip);

a = -(1-E1)*E2*sf;
b = (1-E1)*sf;
c = E2*(E1-1)*(1+cf);
d = 1-E1*cf-(E1-cf)*E2^2;

s0 = (a*exp(i*phi)+b) / (c*cos(phi)+d);

s = s0*exp(-TE/T2)*exp(-i*phi*(TE/TR));
 






