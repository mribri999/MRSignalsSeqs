%	function im = ft(dat)
%
%	Function does a centered 2DFT of the data:
%
%	im = fftshift(fft2(fftshift(dat)));

function im = ft(dat)

im = fftshift(fft2(fftshift(dat)));

