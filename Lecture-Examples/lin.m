%  This function linearizes a matrix (m) of any dimension (eg M=m(:)).
%  If an index vector (ind) is given then the the ind entries of m(:) are
%  returned.
%
% SYNTAX: m=lin(m);
%         m=lin(m,ind);
%
% DBE 12/22/03

function m=lin(m,ind);

m=m(:);

if nargin==2
  m=m(ind);
end

return