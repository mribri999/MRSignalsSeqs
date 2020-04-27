% function im = diamond(size,w)
%
% Make a diamond of given width in a given image size
%
%	INPUT:
%		size = matrix size (2D)
%		w = width (widest point)
%
%	OUTPUT:
%		image that is 1 inside diamond, 0 elsewhere.
%

function im = diamond(size,w)

[x,y] = meshgrid([1:size]-size/2,[1:size]-size/2);
im = ones(size,size);
im(find(x+y>w/2))=0;
im(find(x-y>w/2))=0;
im(find(-x+y>w/2))=0;
im(find(-x-y>w/2))=0;




