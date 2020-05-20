% 	function im = psf2d(ksp,oversamp,timemode)
%
%	Function calculates and plots the point-spread function for
% 	a given 1D or 2D k-space.
%
%	INPUT:
%		ksp = k-space sample points (can be weighted if desired)
%		oversamp = oversampling of PSF (to smooth).  Default 16
%		timemode = 1 for ksp to be k x time
%
%	OUTPUT:
%		im = 2D PSF (can be plotted)
%		plots of PSF.

function im = psf2d(ksp,oversamp,timemode)

if (nargin < 1) ksp = ones(256,256); end;
if (nargin < 2) oversamp=16; end;
if (nargin < 3) timemode=0; end;

kplot=ksp;
if (timemode==1)
  ksp = ifftshift(ifft(ifftshift(ksp,2),[],2),2);	% FFT to interpolate
end;

[y,x] = size(ksp);
ki = ksp;
ki(oversamp*y,oversamp*x)=0;
ki = circshift(ki,round((oversamp/2-0.5)*[y x]));

im = ft(ki);			% 2DFT
im = im/max(abs(im(:)));	% Normalize

subplot(2,2,1);
dispim(kplot); title('k-space'); axis off;

subplot(2,2,2);
dispim(cropim(kplot,16)); title('Central k-space'); axis off;

subplot(2,2,3);
dispim(log(1+abs(im))); title('PSF Image (log scale)'); axis off;

subplot(2,2,4);
npix=8; tt = sprintf('Central PSF %d x %d pixels',npix,npix);
cim = cropim(im,npix*oversamp);
dispim(cim); title(tt); hold on; axis off;
h=plot([1:npix*oversamp],npix*oversamp/2*(1-abs(cim(4*oversamp+1,:))),'c--');
set(h,'LineWidth',2');
h=plot(npix*oversamp/2*(1+abs(cim(:,4*oversamp+1))),[1:npix*oversamp],'y--');
set(h,'LineWidth',2');
% -- Plot 'axis'
plot([1:npix*oversamp],npix/2*oversamp*ones(npix*oversamp),'k-');
plot([1:npix*oversamp],npix/2*oversamp*ones(npix*oversamp),'w--');
plot(npix/2*oversamp*ones(npix*oversamp), [1:npix*oversamp], 'k-');
plot(npix/2*oversamp*ones(npix*oversamp), [1:npix*oversamp], 'w--');
hold off;









