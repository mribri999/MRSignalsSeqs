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
subplot(3,3,1);
plot(lpf); lplot('ky','LPF','a) LPF Magnitude');
subplot(3,3,2);
plot(mf); lplot('ky','Merging Filter','b) Merge Filter');
setprops;

% Low-res Phase Estimate
%
ksp = ifftshift(ifft2(ifftshift(im)));
klp = diag(lpf)*ksp;
imlp= ft(klp);					% Low-Res image
subplot(3,3,4);
dispim(imlp); axis off; title('c) Low-Res Image');
subplot(3,3,5);
dispangle(imlp); axis off; title('d) Phase');

%% f) Phase estimate
khf = diag(mf)*ksp;			% Half kspace
imhf = ft(khf);				% Half Fourier Image (zero filled)
subplot(3,3,7);
dispim(imhf); axis off; title('e) Zero-Filled Image');
subplot(3,3,8);
dispangle(imhf); axis off; title('f) Phase');

phest = angle(imlp);
imhd = imhf.*exp(-i*phest);	% Homodyne recon Image

imhd = imhd * (mean(abs(im(:)))/mean(abs(real(imhd(:))))); % Normalize

subplot(3,3,3);
[l,h] = dispim(real(im)); axis off; title(' Original Image');

subplot(3,3,6);
dispim(real(imhd)); axis off; title(' Final Image');

subplot(3,3,9);
dispim(real(imhd)-abs(im),l,h); axis off; title('Difference Image');


