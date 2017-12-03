% Function imhd = homodyneshow(im,w)
%
%	Function does a homodyne reconstruction on
%	data for a complex image im, with merge filter
%	width w (even) and displays stages.
%

function imhd = homodyneshow(im,w)

if (nargin < 2) w=16; end;
% Get image size
[m,n] = size(im);

% Low-pass filter
%
lpf = zeros(m,1);
lpf( floor(m/2)-w/2: floor(m/2)+w/2-1 ) = 1;	% LPF of 1s in center.

% Merging filter
mf = cumsum(lpf);				% Ramp merge filter
mf = 2*mf/max(mf(:));				% Scale from 0 to 2.
subplot(3,2,1);
plot(lpf); lplot('ky','LPF','LPF Magnitude');
subplot(3,2,2);
plot(mf); lplot('ky','Merging Filter','Merge Filter');
setprops;

% Low-res Phase Estimate
%
ksp = ifftshift(ifft2(ifftshift(im)));
klp = diag(lpf)*ksp;
imlp= ft(klp);					% Low-Res image
subplot(3,2,3);
dispim(imlp); axis off; title('e) Low-Res Image');
subplot(3,2,4);
dispangle(imlp); axis off; title('Phase');

%% f) Phase estimate
khf = diag(mf)*ksp;			% Half kspace
imhf = ft(khf);				% Half Fourier Image (zero filled)
subplot(3,2,5);
dispim(imhf); axis off; title('f) Zero-Filled Image');
subplot(3,2,6);
dispangle(imhf); axis off; title('Phase');

phest = angle(imlp);
imhd = imhf.*exp(-i*phest);	% Homodyne recon Image


