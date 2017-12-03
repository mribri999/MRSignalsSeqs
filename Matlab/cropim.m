%	Function [cim] = cropim(im,sx,sy)
%
%	Crops image to sx x sy.
%		If sy is omitted, sy=sx.
%		If sx is omitted, sx=size(1)/2.
%		If sx or sy is < 1, then fractional
%

function [cim] = cropim(im,sx,sy)

sz = size(im);
if (nargin < 2) sx = floor(sz(1)/2); end;
if (nargin < 3) sy = sx; end;

if (sx<1) sx=sz(1)*sx; end;	% If fractional
if (sy<1) sy=sz(1)*sy; end;


stx = floor(sz(2)/2-sx/2)+1;
sty = floor(sz(1)/2-sy/2)+1;

cim = im( sty:sty+sy-1, stx:stx+sx-1 );


