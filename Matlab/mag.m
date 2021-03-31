% This function finds the magnitude of a vector of any length or of a
% matrix along a given dimension (1 or 2)
%
% SYNTAX: mag=mag(vector,dim);
%         mag=mag(vector);     % If a vector 
%
% DBE 01/25/00

function mag=mag(vector,dim)

% len=length(vector);
if nargin==1
  dim=find(min(size(vector))==size(vector));
end

if     size(vector,1)==1
  dim=2;
elseif size(vector,2)==1
  dim=1;
end

mag= sqrt(sum(vector.^2,dim));

return;