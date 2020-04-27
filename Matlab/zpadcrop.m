%	function newim = zpadcrop(im,newsize)
%
%	Function zero pads or crops an image to the new
%	size in x and y, keeping the image centered.
%

function newim = zpadcrop(im,newsize)


sz = [size(im)];
if (length(sz)==2) sz = [sz 1]; end;	% Make 3D at least.

 
newim = zeros(newsize);


if (length(newsize) < length(sz))
	% Copy extra dimensions from original
	newsize(length(newsize)+1:length(sz))= sz(length(newsize)+1:end);
else
	% Crop extra dimensions in new size
	newsize = newsize(1:length(sz));
end;

cz = floor(sz/2)+1;		% Central point.
ncz= floor(newsize/2)+1;	% New central point

minsz = min([sz; newsize]);	% minimum of sizes

% -- Start and end indices
newst = ncz-floor(minsz/2);
newen = ncz+ceil(minsz/2)-1;
oldst = cz-floor(minsz/2);
olden = cz+ceil(minsz/2)-1;


newim(newst(1):newen(1),newst(2):newen(2),newst(3):newen(3)) = im(oldst(1):olden(1),oldst(2):olden(2),oldst(3):olden(3));


