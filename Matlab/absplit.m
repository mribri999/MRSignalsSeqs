
%	function [Ai,Bi] = absplit(n,A,B)
%
%	Function 'splits' a propagation A,B into sub-matrices so that
% 	n 'propagations' of M'=Ai*M+Bi is the same as M'=A*M+B.
%
%	This is useful for animating an operation, where you want
%	to split it into many frames.  See abanim.m
%
function [Ai,Bi] = absplit(n,A,B)

if (nargin < 3) B=[0;0;0]; end;

Ai = A^(1/n);	% Multiplications are just nth root

Aprop = inv(eye(3)-Ai)* (eye(3)-Ai^n);	% Series formula for sum A^k
Bi = inv(Aprop)*B;


