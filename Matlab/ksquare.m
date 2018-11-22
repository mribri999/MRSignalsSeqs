% function kdata = ksquare(center,swidth,kloc,tsamp,df)
%
%	Function generates k-space data for a square 
%	with given center location (cm) and swidth length (cm)
%	at the given k-space locations, with the given sampling period
%	and off-resonance (Hz)
%
%	INPUT:
%		center,swidth = describes square, in cm
%			center can be a list for more squares
%		kloc = array of kx+i*ky, Nsamples x Ninterleaves (cm^-1)
%		tsamp = time between points (resets each interleaf) (sec)
%		df = off-resonance (Hz)
%	OUTPUT:
%		Data samples, same size as kloc
%		
function kdata = ksquare(center,swidth,kloc,tsamp,df)

if (nargin < 5) df = 0; end;		% on resonance
if (nargin < 4) tsamp = 0.000004; end;	% 4us
if (nargin < 3) 
  [kx,ky] = meshgrid([-128:127]/128*5,[-128:127]/128*5);
  kloc = kx+i*ky;
end;
if (nargin < 2) swidth = 1.9; end;
if (nargin < 1) center = 0; end;

% -- Define data with a sinc
sdata = swidth*sinc(swidth*real(kloc)).*sinc(swidth*imag(kloc));
kdata = 0*sdata;

% -- Add linear phase for shift of center
for q=1:length(center)
  
  thisk = exp(i*2*pi*real(center(q))*real(kloc)).*sdata;
  thisk = exp(i*2*pi*imag(center(q))*imag(kloc)).*thisk;
  kdata = kdata + thisk;
end;

% -- Add linear phase for off-resonance.
if (nargin > 4)
  ph = exp(2*i*pi*tsamp*[1:size(kloc,1)]*df);	% phase due to off-resonance.
  kdata = diag(ph)*kdata;			% off-res induced phase
end;







