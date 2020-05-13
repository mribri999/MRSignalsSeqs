% Lecture 11, Example 01
%
% Circular vs square PSF

N = 256;	% Image
No= 32;		% oversample to interpolate PSF

kcirc = circ(N*No,N/2);
krect = 0*kcirc;
krect(N*No/2-N/2+1:N*No/2+N/2, N*No/2-N/2+1:N*No/2+N/2 )=1;

% -- FFTs
prect = ft(krect);
pcirc = ft(kcirc);

% -- Normalize & Crop
prect = cropim(prect / max(abs(prect(:))),4*No,4*No);
pcirc = cropim(pcirc / max(abs(pcirc(:))),4*No,4*No);

set(0,'defaultAxesFontSize',18);        % Default font sizes)

subplot(1,3,1);
dispim(prect);
dispim(prect,0,1); title('Square Sampling PSF');
axis off;

subplot(1,3,2);
dispim(pcirc,0,1); title('Circular Sampling PSF');
axis off;

subplot(1,3,3);
dispim(prect-pcirc,0,1); title('Difference in PSFs');
axis off;
colorbar;

colormap('jet');

