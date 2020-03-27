function plotm(M)
%function plotm(M)
%	Show spins on a plot with 3D, Mx-My, Mx-Mz, My-Mz views.
%
%	INPUT:
%		M = 3xN magnetization to plot
%		scale = (optional) scale for plots
%
%	OUTPUT:
%		plot on current figure.


if (nargin < 1) 
  disp('No argmenents passed to plotm - using defaults')	
  M = [0.8,0,0.2].'; 
end;

clf;
subplot(2,2,1);
showspins(M);
title('3D View')

subplot(2,2,2);
showspins(M);
view(0,90);
tt = sprintf('Mx - My'); title(tt);

subplot(2,2,3);
showspins(M);
view(0,0);
tt = sprintf('Mx - Mz'); title(tt);

subplot(2,2,4);
showspins(M);
view(90,0);
tt = sprintf('My - Mz'); title(tt);

setprops;

