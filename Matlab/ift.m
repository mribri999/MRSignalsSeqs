%	function im = ft(dat)
%
%	Function does a centered inverse 2DFT of the data:
%
%	im = ifftshift(ifft2(ifftshift(dat)));

function im = ift(dat)

im = ifftshift(ifft2(ifftshift(dat)));

