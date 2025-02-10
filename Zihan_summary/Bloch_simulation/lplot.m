%	Function lplot(xlab,ylab,tit,ax,grid);
%
%	Label a figure, turn on grid, and set properties
%

function lplot(xlab,ylab,tit,ax)

if (nargin >0) && (length(xlab)>0)
  xlabel(xlab); end;

if (nargin >1) && (length(ylab)>0)
  ylabel(ylab); end;

if (nargin >2) && (length(tit)>0)
  title(tit); end;

if (nargin >3) && (length(ax)>3)
  axis(ax); end;

if (nargin <5) || (grid==1)
  grid on; end;

setprops;

