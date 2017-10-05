% Calculation of simple 1D diffusion sensitivity.

D = 1e-6;	% mm^2/ms
G = 0.040;		% mT/mm;
T = 10;			% ms.
gamma = 42.58	% kHz/mT
dt = 20;	% ms


sig = sqrt(2*D*dt);	% mm.

x = [-3:.01:3]*sig;	% positions.

% == Numerical Integration to test 
sum = 0;
dx = x(2)-x(1);
for k=1:length(x);
  val = cos(2*pi*gamma*G*T*x(k)) * 1/sqrt(4*pi*D*dt) * exp(-x(k)^2/(4*D*dt));
  sum = sum+val*dx;
end;

sum

% == Calculated Solution
sig = exp(-(2*pi*gamma*G*T)^2*dt * D)
