%	Short-TR IR Signal Example
%
%	Use abprop to do this compactly.
%	Note we start applying A,B from TE onward as this is steady-state.
%

% Simple Inversion-Recovery sequence where TR is not long enough for
% full recovery so steady-state must be calculated.

% Times are all in seconds here.
TR=1;
TI = 0.5;
TE = 0.05;
T1 = 0.5;
T2 = 0.1;

[A,B,Mss] = abprop(relax(TR-TE-TI,T1,T2,1),xrot(180), ...
	relax(TI,T1,T2,1),xrot(90), relax(TE,T1,T2,1));

Mss	
