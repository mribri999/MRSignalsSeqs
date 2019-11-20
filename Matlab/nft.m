%	function im = nft(dat)
%
%	Function does a centered, normalized 2DFT of the data:
%
%	im = fftshift(fft2(fftshift(dat)))/sqrt(length(dat(:)));

function im = nft(dat)

im = fftshift(fft2(fftshift(dat))) / sqrt(length(dat(:)));

