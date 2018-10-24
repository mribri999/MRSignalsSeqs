%	function y = sinc(x)
%
%	Function returns y = sin(pi*x)/pi*x if x is non-zero, and 1 otherwise.
%	Works on arrays.

function y = sinc(x)


sz = size(x);

x = x(:);
y = x;			% Allocate
f1 = find(x==0);
f2 = find(x~=0);

y(f1)=1;
y(f2) = sin(pi*x(f2))./(pi*x(f2));

y = reshape(y,sz);



