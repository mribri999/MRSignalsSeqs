function plotm(M,scale,orig)
%function plotm(M,scale,orig)
%	Show spins on a plot with 3D, Mx-My, Mx-Mz, My-Mz views.
%
%	INPUT:
%		M = 3xN magnetization to plot
%		scale = (optional) scale for plots
%		orig = 3xN origins for vectors
%
%	OUTPUT:
%		plot on current figure.


if (nargin < 1) 
  disp('No argmenents passed to plotm - using defaults')	
  M = [0.8,0,0.2].'; 
end;
if (nargin < 2) scale = 1; end;
if (nargin < 3) orig = 0*M; end;

M=real(M);	% Avoid errors! 
clf;
subplot(2,2,1);
showspins(M,scale,orig);
title('3D View')

subplot(2,2,2);
showspins(M,scale,orig);
view(0,90);
tt = sprintf('Mx - My'); title(tt);

subplot(2,2,3);
showspins(M,scale,orig);
view(0,0);
tt = sprintf('Mx - Mz'); title(tt);

subplot(2,2,4);
showspins(M,scale,orig);
view(90,0);
tt = sprintf('My - Mz'); title(tt);

setprops;

