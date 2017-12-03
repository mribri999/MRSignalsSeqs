%
%	Function [dat] = gridmat(ksp,kdat,dcf,gridsize)
%
%	function performs gridding in Matlab.  This is designed to
%	be a teaching function, and is not necessarily fast(!), but
%	reasonable for 2D gridding.
%
%	INPUT:
%         ksp = k-space locations, kx+i*ky, normalized to |k|<=0.5, 1/pixel.
%         dat = complex data samples at k-space locations.
%	  dcf = density compensation factors at k-space locations.
%	  gridsize = size of grid (default is 256)
%
%	OUTPUT:
%	  dat = matrix of gridded data samples
%
%	Many parameters come from Jackson 1991 here.

%	B.Hargreaves - October 2014.

function [dat] = gridmat(ksp,kdat,dcf,gridsize)

ksp=ksp(:);		% Translate to 1D vector in case this is 2D
dcf=dcf(:);
kdat=kdat(:);

kwid = 3;		% Convolution kernel width
kbbeta = 4.2;		% Kaiser-bessel Beta

% -- Allocate grid to accumulate data
padgridsize = gridsize+4*kwid;			% Padded grid to avoid errors
padgrid = zeros(padgridsize,padgridsize);	% Padded grid.

% -- Sample Density correction
kdat = kdat .* dcf;		% Density correct (Simple!).

% -- Use a lookup-table for grid kernel
dr = 0.01;			 
r = [0:dr:2*kwid];		% input to kaiser-bessel
kerntab = kb(r,kwid,kbbeta);	% kaiser-bessel function

% -- Scale kspace points to grid units
kxi = real(ksp)*(gridsize-1)+padgridsize/2;	% Scale kx to grid units
kyi = imag(ksp)*(gridsize-1)+padgridsize/2;	% Scale ky to grid units

% -- Round to nearest grid point
ikxi = round(kxi);	% Closest integer
ikyi = round(kyi);	% Closest integer

% -- Make a small matrix, that completely encompasses the grid points
% 	within the convolution kernel around a k-space sample --

sgridext = ceil(kwid/2)+1;	% Size of small grid around sample
[smaty,smatx] = meshgrid([-sgridext :sgridext ],[-sgridext :sgridext ]);
sgrid = 0*smatx;		% allocate


% -- Go through k-space samples to do convolution

for p = 1:length(ksp)
  %tt = sprintf('Gridding sample %d of %d',p,length(ksp));
  %disp(tt);

  gridx = smatx + ikxi(p);	% grid of 'closest' integer pts
  gridy = smaty + ikyi(p);	% same in y	
			
  % -- Distance index (array), to use for kernel lookup	 
  %    Just calculating kernel value between ksp(p) and every point
  %    in this small grid.
  dist = round(sqrt((gridx - kxi(p)).^2 + (gridy -kyi(p)).^2)/dr)+1;

  sgrid(:) = kdat(p) * kerntab(dist(:)); 	% Convolve sample w/ kernel

  % -- Add the 'small-grid' into the padded grid
  padgrid(ikxi(p)-sgridext:ikxi(p)+sgridext,ikyi(p)- ...
	sgridext:ikyi(p)+sgridext) = padgrid(ikxi(p)-sgridext: ...
	ikxi(p)+sgridext,ikyi(p)-sgridext:ikyi(p)+sgridext) + sgrid;
   
end;

% -- Extract the main grid from the padded grid.
dat = padgrid(2*kwid+1:2*kwid+gridsize, 2*kwid+1:2*kwid+gridsize);






% Kaiser-bessel Kernel


%
%	function y = kb(u,w,beta)
%
%	Computes the Kaiser-Bessel function used for gridding, namely
%
%	y = f(u,w,beta) = I0 [ beta*sqrt(1-(2u/w)^2) ]/w
%
%	where I0 is the zero-order modified Bessel function
%		of the first kind.
%
%	INPUT:
%		u = vector of k-space locations for calculation.
%		w = width parameter - see Jackson et al.
%		beta = beta parameter - see Jackson et al.
%
%	OUTPUT:
%		y = vector of Kaiser-Bessel values.
%
%	SEE ALSO:
%		kbk2x.m

%	B. Hargreaves	Oct, 2003.




function y = kb(u,w,beta)

if (nargin < 3)
	error('Not enough arguments -- 3 arguments required. ');
end;


if (length(w) > 1)
	error('w should be a single scalar value.');
end;


y = 0*u;				% Allocate space.
uz = find(abs(u)< w/2);			% Indices where u<w/2.

if (length(uz) > 0)				% Calculate y at indices uz.
	x = beta*sqrt(1-(2*u(uz)/w).^2);	% Argument - see Jackson '91.
	y(uz) = besseli(0,x)./w;
end;

y = real(y);		% Force to be real.




