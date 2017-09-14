%function [FpFmZ] = epg_spins2FZ(M,trim)
%
%	Function converts a 3xQ array of vectors to EPG states F+,F- and Z,
%	with Q representing the dimension across a voxel.  Note that it
%	is not realistic to have more than floor(Q/2)+1 states, since M
%	is real-valued and the F states can be complex.
%
%	INPUT:
%		M = [Mx My Mz].' representing magnetization vectors.
%		trim = threshold to trim, 0 for none (default = 0.01)
%		
%	OUTPUT:
%		FpFmZ = F,Z states, 3xQ vector where there are Q states,
%			where rows represent F+, F- and Z states starting at 0.
%
function [FpFmZ] = epg_spins2FZ(M,trim)

if (size(M,1)~=3) error('EPG state must be 3xQ vector'); end;
if (nargin < 2) trim=0.01; end;
Q = size(M,2);

% -- The following are Eqs. 2-4 from Wiegel 2010:
Fp = fft(M(1,:)+i*M(2,:));		% FFT to get F+ states.
Fm = fft(M(1,:)-i*M(2,:));		% FFT to get F- states.
Z = fft(M(3,:));			% FFT to get Z.

FpFmZ = [Fp;Fm;Z];			% Combine into 3xQ matrix
FpFmZ = FpFmZ(:,1:floor(Q/2)+1);	% Cut off redundant states.
FpFmZ = epg_trim([FpFmZ],trim);		% Trim near-zero states.



