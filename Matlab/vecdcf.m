%	function [dcf,k] = vecdcf(g,k [,N,res,FOV])
%
%	Calculate the Density Correction Factors for
%	a given spiral trajectory using the vector approach
%	used by Craig Meyer (MRM, 1992)
%
%	INPUT:
%		g = gradient (Gx + i*Gy)
%		k = k-space  (kx + i*ky)
%
%
%	OUTPUT:
%		dcf = list of density compensation factors
%
%
%	NOTES:
%		if k is omitted, it is calculated from g.
%		If k is 1x1, it is assumed to be the
%			#of sample points.
%		

% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: vecdcf.m,v $
%	Revision 1.1  2003/08/28 20:27:01  brian
%	Added.
%	
%	Revision 1.1  2002/04/25 21:58:20  bah
%	Added to CVS.
%	
%
% ===========================================================




function [dcf,k] = vecdcf(g,k,N,res,FOV)


% ========== Fix Gradient if Nx2 ============
sg = size(g);
if (sg(2)==2)		% Gradient vector is Nx2
	g = g(:,1)+i*g(:,2);
end;

% ========== Calculate k-space from gradients ============
if ((nargin < 2) | (length(k)==0))	% If only g is given.
	k = cumint(g,0.5);
	k = k*0.5/max(abs(k(:)));
end;



if (max(size(k))==1)
	Nsamps = k;
	k = cumint(g(1:Nsamps),0.5);
	k = k*0.5/max(abs(k(:)));
end;

Nsamps = length(k);

	

% ========== Calculate Density Corr Factors using Vector Approach =====

g3 = [real(g(:)) imag(g(:)) 0*g(:)]';
k3 = [real(k(:)) imag(k(:)) 0*k(:)]';
dcf = 0*k;
for m=1:Nsamps
	dcf(m) = norm(cross(g3(:,m),k3(:,m)));
	if (abs(k(m)) > 0)
		dcf(m) = dcf(m)/abs(k(m));
	else
		dcf(m) = 0;
	end;
end;


% =========== IF N, res, FOV given, replicate for N interleaves ===
dcf = dcf(:);
k = k(:);

if (nargin > 2)
	dcfall = zeros(length(dcf),N);
	kall = zeros(length(dcf),N);
	for p=1:N
		ph=exp(i*2*pi*(p-1)/N);
		kall(:,p)=k*ph;
		dcfall(:,p)=dcf;
	end;
dcf=dcfall;
k = kall * FOV/res/256;
end;





