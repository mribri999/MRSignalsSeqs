% [sig] = gresignal(T1,T2,TE,TR,flip)
%
%	Plot the theoretical signal for a gradient-spoiled sequence.
%	Note that this is NOT RF-spoiled GRE, which is the same
%	(roughly) as excitation-recovery.
%
%	T1,T2,TE,TR are parameters for tissue and sequence.
% 	flip is the flip angle in degrees.
% The signal as a fraction of Mo is returned.
%
%	See Buxton 1989.


function [sig] = gresignal(T1,T2,TE,TR,flip)

E1 = exp(-TR/T1);
E2 = exp(-TR/T2);
E  = exp(-TE/T2);

s = sin(pi*flip/180);
c = cos(pi*flip/180);

B = 1 - E1*c - E2^2*(E1-c);
C = E2*(1-E1)*(1+c);
A = sqrt(B^2-C^2);

if (abs(A*C)==0)
  sig = 1;			
else
  sig = s*(1-E1)*E* (C-E2*(B-A)) / (A*C);
end;


