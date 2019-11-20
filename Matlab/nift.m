%	function im = nift(dat)
%
%	Function does a centered, normalized inverse 2DFT of the data:
%
%	im = ifftshift(ifft2(ifftshift(dat))) * sqrt(length(dat(:)));

function im = nft(dat)

im = ifftshift(ifft2(ifftshift(dat))) * sqrt(length(dat(:)));

