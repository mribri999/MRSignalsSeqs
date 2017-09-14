%
%	function [g,t] = mingrad(area,Gmax,Smax,dt)
%
%	Function calculates the fastest gradient waveform
%	to achieve the given area in k-space, and returns
%	it, along with a time array.
%
%	INPUT:
%		area = desired k-space area in cm^(-1)
%		Gmax = maximum gradient (mT/m) [50];
%		Smax = maximum slew (mT/m/ms)  [200];
%		dt = time step in ms [.004];
%
%	OUTPUT:
%		g = gradient waveform (mT/m);
%		t = time array (ms);
%
function [g,t] = mingrad(area,Gmax,Smax,dt)

area = area*100;			% Convert cm(-1) to m(-1)
if (nargin < 2) Gmax = 50; end;		% mT/m
if (nargin < 3) Smax = 200; end;	% mT/m/ms
if (nargin < 4) dt = 0.004; end;	% ms

% 2 cases - triangle or trapezoid.
%
gamma = 42.58;	% kHz/mT
gamma = 40;	% kHz/mT

% Area of largest triangle
Atri = gamma * (Gmax^2/Smax); 	% 1/2 base * height


if (area <= Atri)			% Area = (tramp)^2*Smax
  tramp = sqrt(area/gamma/Smax); 		% Ramp time (ms)
  Nramp = ceil(tramp/dt);  			% #pts on ramp
  g = [1:Nramp]*Smax*dt;			% Generate ramp
  g = [g fliplr(g)];
else
  tramp = (Gmax/Smax);				% Ramp time (ms)
  Nramp = ceil(tramp/dt);  			% #pts on ramp
  gramp = [1:Nramp]/Nramp*Gmax;
  Nplat = ceil(area/gamma/Gmax/dt - Gmax/Smax/dt);
  g = [gramp Gmax*ones(1,Nplat) fliplr(gramp)];
end;

t = [1:length(g)]*dt;			% Generate time array
g = g * (area/gamma)/(sum(g)*dt);	% Correct for rounding.



