function [p,arrc] = aclinphase(arr,pwr)
%function [p,arrc] = ahncholinphase(arr,pwr)
%
%       Calculate the "Ahn-Cho" linear phase across 
%	a complex array, per pt.
%
%	arr = complex data
%	pwr = power to which to raise input to make more noise immune.
%		A very small pwr will weight points more equally (ex 0.1)
%
%	p = linear phase in radians per pixel.

if (nargin < 2) pwr=2; end;
a1 = arr(1:end-1,:,:,:,:).^pwr; % Raise to pwr to make more noise-immune.
a2 = arr(2:end,:,:,:,:).^pwr;

p = angle(sum((a2.*conj(a1)),1))/pwr;

if (nargout >1)
  arrc = arr(:) .* exp(-i*p*[-floor(length(arr)/2):floor((length(arr)-1)/2)].');
  arrc = reshape(arrc,size(arr));
end;

