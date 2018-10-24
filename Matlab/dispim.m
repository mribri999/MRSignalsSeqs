% dispim(im,[low,high])
% displays magnitude version of image im (complex array)
%
%	im = 2D image array (real or complex, magnitude displayed)
%	low = black level, high = white level.
%	Omit low,high for autoscaling.
%
% ===========================================================


function [low,high] = dispim(im,low,high)

im = squeeze((im));

if (nargin < 2)
	low = 0;	
end;
if (nargin < 3)
	immax = max(abs(im(:)));
	imstd = std(abs(im(:)));
	high = immax-.5*imstd;
end;


% display image
scale= 256/(high-low);
offset = scale*low;

% set colormap for 256 gray levels
c = colormap;
if (sum(abs(std(c')))>0.01)	% Only change colormap if it is not gray.
	a=[1 1 1]/256;
	b=[1:256];
	c=(a'*b)';
	colormap(c);
end;

image(abs(im)*scale-offset);
axis('square');


