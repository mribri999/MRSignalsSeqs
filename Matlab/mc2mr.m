%
%	function Mr = mc2mr(Mc)
%
%	Transform magnetization of the 
%	complex form [mx+i*my; mx-i*my; mz].
%	to the rectilinear  form [mx;my;mz]. to 
%	INPUT:
%		3xN array of rectilinear magnetization [mx;my;mz]
%
%	OUTPUT:
%		3xN array of complex magnetization [mx+i*my; mx-i*my; mz]
%
function Mr = mc2mr(Mc)

% -- Check input is 3xN, add rows if not.
sz = size(Mc);
if (sz(1)~=3) error('Mc should have 3 rows'); end;

T = [1,1,0;-i,i,0;0,0,2]/2;	% Transformation
Mr = T*Mc;

% -- Alternative method 
%Mr = 0*Mc;			% Allocate
%Mr(1,:) = real((Mc(1,:)+Mc(2,:))/2);  % 1st row
%Mr(2,:) = imag((Mc(1,:)-Mc(2,:))/2);  % 2nd row
%Mr(3,:) = Mc(3,:);


