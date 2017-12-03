%
%	Function ph = lfphase(m,n,w)
%
%	Function generates a low-resolution phase function, by
%	FFT of random coefficients.
%
%	2*w+1 = width of non-zero Fourier coefficients.

function ph = lfphase(m,n,w)

if (nargin < 2) n=m; end;
if (nargin < 3) w=4; end;

arr = zeros(m,n);
arr(floor(m/2-w):floor(m/2+w) ,floor(n/2-w):floor(n/2+w)) = randn(2*w+1,2*w+1);
ang=fftshift(fft2(fftshift(arr)));     % Angle
ang= real(ang);
ang= 2*pi*(ang-min(ang(:))) / (max(ang(:))-min(ang(:))) -pi;
ph = exp(i*ang);



