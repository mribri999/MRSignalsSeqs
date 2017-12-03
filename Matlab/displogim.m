% displogim(im)
% displays log-magnitude version of image im (complex array)
%
%	im = 2D image array (real or complex, magnitude displayed)
%
% ===========================================================


function displogim(im)

im = squeeze((im));
lowexp = 0;		% Make negative if max(im) is small.

im = log(abs(im));
f = find(im(:)< (lowexp));
im(f)= lowexp;

im = im - lowexp;
dispim(im);

