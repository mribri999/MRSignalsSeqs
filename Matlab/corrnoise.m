%	Function [nc,Rout] = corrnoise(mn,R,n)
%
%	Function generates correlated noise, where
%	R can be complex.
%
%	INPUT:
%		mn = mean of noise (default to 0)
%		R = noise covariance matrix
%		n = number of samples
%
%	OUTPUT:
%		nc = mxn array of noise samples, where m is dimension of R
%		Rout = covariance of samples (check!)

%	Taken from https://www.mathworks.com/matlabcentral/fileexchange/21156-correlated-gaussian-noise

function [nc,Rout] = corrnoise(mn,R,n)

if (nargin < 3)
  n = 100;
end;
if (nargin < 2) || (size(R,1)<1)
  R = [1 0.1 0.2*i; 0.1 0.8 0.3; -0.2*i 0.3 0.95];
  %R = [1 0.1 0.2; 0.1 0.8 0.3; 0.2 0.3 0.95];
end;
if (nargin < 1) || (length(mn)<1)
  mn = zeros(size(R,1),1);
end;

nvect = randn(size(R,1),n);	% Noise vector

% symmetrize matrix R
R = 0.5 * (R + R');		% Force correlation matrix to be symmetric.

R

[v,d] = eig(R);			% Eigen decompolstion.

if (any(diag(d) <=0))
	error('R must be positive definite');
end;

w = v*sqrt(d);

nc = w*nvect;


if (nargout > 1)
  Rout = (nc * nc')/n;
end;
 






