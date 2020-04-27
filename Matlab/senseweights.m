%
% function [W,gfact,imnoise] = senseweights(sens,Psi)
%
%	Function calculates the SENSE weights at a pixel given
%	input sensitivities and covariance matrix.  A reduction
%	factor R can be assumed from the shape of the sens matrix.
%	The number of coils Nc is determined by the size of sens.  If Psi
%	is omitted it is assumed to be the identity matrix.
%
%	INPUT:
%		sens = Nc x R matrix of coil sensitivites 
%		Psi = Nc x Nc noise covariance matrix.
%
%	OUTPUT:
%		W = coil combination weights - R x Nc
%		gfact = gfactor (1xR vector including each aliased pixel)
%		imnoise = image noise covariance matrix (RxR)
%
function [W,gfact,imnoise] = senseweights(sens,Psi)

[Nc,R] = size(sens);		   % -- Get #coils and reduction factor
if nargin < 2 Psi = eye(Nc); end;  % Default to perfectly uncorrelated noise

% --
invPsi = inv(Psi);		% Calculate inverse once
CiPsiC = sens'*invPsi*sens; 		% Calculate jsut once
iCiPsiC = inv(CiPsiC);		% Calculate just once

W = iCiPsiC * sens'*invPsi;	% SENSE weights
imnoise = iCiPsiC;		% Image noise covariance
gfact = sqrt(diag(imnoise).* diag(inv(imnoise)));	% Elementwise prod.


