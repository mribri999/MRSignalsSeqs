function [arr] = circ(size,r)
%function [arr] = circ(size,r)
%	circ	- generates an array which contains a
%		cirle function, that is 1 for points
%		within r of the center of the array and
%		0 elsewhere.
%
%		size = array dimension.
%		r = circle radius.
%

% 	Brian Hargreaves, Aug 1/96.
%

arr = zeros(size,size);
x = 1:size;
c = (size+1)/2;
x = x-c;

arrx = ones(size,1)*x;
arry = arrx';
r2 = r^2;
arrr = sqrt(arrx.*arrx + arry.*arry);
f = find(arrr(:)<r);
arr(f) = 1;


    


