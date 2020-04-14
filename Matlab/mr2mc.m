%
%	function Mc = mr2mc(Mr)
%
%	Transform magnetization of the form [mx;my;mz] to 
%	complex form [mx+i*my; mx-i*my; mz].
%
%	INPUT:
%		3xN array of rectilinear magnetization [mx;my;mz]
%		(if 1xN or 2xN, will zero-pad)
%
%	OUTPUT:
%		3xN array of complex magnetization [mx+i*my; mx-i*my; mz]
%
function Mc = mr2mc(Mr)

% -- Check input is 3xN, add rows if not.
sz = size(Mr);
if (sz(1)~=3) Mr(3,1)=0; Mr=Mr(1:3,:); 
disp('mr2mc.m:  Warning input is not 3 rows - zero-filling and cropping'); end;  
T = [1,i,0;1,-i,0;0,0,1];     % Transformation Mr to Mc
Mc = T*Mr;

% -- Alternative method
Mc = 0*Mr;			% Allocate
Mc(1,:) = Mr(1,:)+i*Mr(2,:);	% 1st row
Mc(2,:) = Mr(1,:)-i*Mr(2,:);	% 2nd row
Mc(3,:) = Mr(3,:);		% mz is the same.

