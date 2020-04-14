function FZ = epg_m0(N)
% function FZ = epg_m0(N)
%
%	Function returns an EPG state matrix with order N-1
%
if (nargin < 1) N=1; end;
FZ = zeros(3,N);
FZ(3,1)=1;


