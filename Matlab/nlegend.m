%	Function nlegend(a,fmt);
%
%	Make a legend from the numbers in an array a.
%	Default is to use %d, but if fmt is defined, use it instead.
%
%	Example if a=[1.2 3.4] and fmt='b=%g' then this will
%	effectively do legend('b=1.2','b=3.4')
%
%	Note that you have to make sure length(a) does not exceed the
%	number of plot traces!
%
function nlegend(a,fmt)

if (nargin < 2) fmt = '%d'; end;

allstr = {};
for k=1:length(a)
  tt = sprintf(fmt,a(k));
  allstr{k} = tt;
end;

legend(allstr);


