%
%	function [A,B] = relax(T,T1,T2,combine)
%
%	Function returns A matrix and B vector for
%	relaxation over time T, with relaxation times T1 and T2,
%	such that m' = A*m+B where m and m' are the magnetization before
%	and after the relaxation.
%
%	If combine==1 and one output, then output is [A B] (3x4)
%
%
function [A,B] = relax(T,T1,T2,combine)

E1 = exp(-T/T1);
E2 = exp(-T/T2);

A = diag([E2 E2 E1]);
B = [0;0;1-E1];

if ((nargin > 3) && (nargout < 2) && (combine==1))
  A = [A B];
end;

  
