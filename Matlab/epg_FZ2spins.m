%function [M] = epg_FZ2spins(FpFmZ, N, frac)
%
%	Function returns a 3xN array of magnetization vectors [Mx My Mz].'
%	that is represented by the EPG state FpFmZ.  Note that the 
%	magnetization assumes that a "voxel" is divided into N states,
%	as is necessary to make the conversion.  N must be at least 
%	as large as twice the number of F+ states minus one.
%
%	INPUT:
%		FpFmZ = F,Z states, 3xQ vector where there are Q states,
%			where rows represent F+, F- and Z states starting at 0.
%		N = number of spins used to represent the states.  Default
%			is minimum.
%		frac = (optional) fraction of state to increment/decrement
%			so if frac=1, this is equivalent to first passing
%			FpFmZ through epg_grad().  This is mostly just 
%			to make plots.
%		
%	OUTPUT:
%		M = [Mx My Mz].' representing magnetization vectors.
%
function [M] = epg_FZ2spins(FpFmZ,N,frac)

if (nargin < 3) frac=0;	end;	% -- Don't include fractional rotation

if (size(FpFmZ,1)~=3) error('EPG state must be 3xQ vector'); end;

Ns = size(FpFmZ,2);	% Number of EPG states (N should be 2*Q-1 or more)
if (nargin < 2) N = 2*Ns-1; end;	% Minimum value for N.

if (N < 2*Ns-1) warning('Number of spins should be at least double number of states minus 1'); end;

			
% -- Use a matrix for FFT to support arbitrary N.
%  	This is because we are going from a discrete set of 
%	coefficients to a "continuous" distribution of Mx,My,Mz
x = [0:N-1]/N-0.5;	% Fraction of a cycle for spins locations.
ph = exp(i*2*pi*x.'*([-(Ns-1):(Ns-1)]+frac));	% Phasors for Fourier Xform
						% Note: frac = 0 usually!



Fstates = [fliplr(conj(FpFmZ(2,2:end))) FpFmZ(1,:)]; 	% Stretched out
					% Vector of Fp (-n,n)
Mxy = ph * Fstates.';			% 1D Fourier Transform to Mxy.

ph = exp(i*2*pi*x.'*[0:Ns-1]);		% Phasors for Mz transform
FpFmZ(3,1)=FpFmZ(3,1)/2;		% Account for discretization
Mz = 2*real(ph*FpFmZ(3,:).');		% Transform to Mz

M = [real(Mxy) imag(Mxy) Mz].'/N;	% Form output vector.


