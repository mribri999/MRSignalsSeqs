% function [im] = dispkspim(ksp)
% displays magnitude/phase of k-space and image in 2x2 array.
%
%	ksp = k-space data
%	im = image
%
% ===========================================================

function im = dispkspim(ksp)

subplot(2,2,1);
dispim(log(abs(ksp)+1));
title('k-space log magnitude');

subplot(2,2,2);
dispangle(angle(ksp));
title('k-space phase');

im = ifftshift(ifft2(ifftshift(ksp)));

subplot(2,2,3);
dispim(im);
title('Image Magnitude');

subplot(2,2,4);
dispangle(im);
title('Image Phase');


