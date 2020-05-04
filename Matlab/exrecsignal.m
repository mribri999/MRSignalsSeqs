%
%[sig] = exrecsignal(T1,T2,TE,TR,flip)
%
% Calculates the steady-state signal of a 
% simple excitation recovery sequence (or SPGR).
%
% T1,T2,TE,TR are tissue/sequence parameters. 
% flip is the flip angle in degrees.
% The signal as a fraction of Mo is returned.
%


function [sig,M] = exrecsignal(T1,T2,TE,TR,flip); 

s = sin(pi*flip/180);
c = cos(pi*flip/180);

E = exp(-TE/T2);
R = exp(-(TR)/T1);

if (abs(1-R*c)<.0001)
  sig = 1;			
  M=[1;0;0];
else
  sig = (1-R)*E*s/(1-R*c);
  M = [sig;0;(1-R)*E*c/(1-R*c)];
end;

