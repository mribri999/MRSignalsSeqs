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
if (sz(1)~=3) Mc(3,1)=0; Mc=Mc(1:3,:); 
disp('mc2mr.m:  Warning input is not 3 rows - zero-filling and cropping'); end; 

Mr = 0*Mc;			% Allocate
Mr(1,:) = real((Mc(1,:)+Mc(2,:))/2);  % 1st row
Mr(2,:) = imag((Mc(1,:)-Mc(2,:))/2);  % 2nd row
Mr(3,:) = Mc(3,:);


