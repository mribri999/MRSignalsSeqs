% function im = gaussian2d(size,mean,sigma)
%
% Function generates a gaussian image with the above parameters.
% Mean and sigma can be vectors if different in x and y, or it
% is the same.

function im = gaussian2d(sz,mn,sig);


if (length(sz) < 2) sz = [sz sz]; end;
im = zeros(sz(2),sz(1));

if (length(mn)<2) mn = [mn mn]; end;
if (length(sig)<2) sig = [sig sig]; end;

x = 1:sz(1);
gx = gaussian(x,mn(1),sig(1));

y = [1:sz(2)]';
gy = gaussian(y,mn(2),sig(2));

im = gy*gx;


