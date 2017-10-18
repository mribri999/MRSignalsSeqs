function epg_show(FZ,frac,scale,Nspins)
%function epg_show(FZ,frac,scale,Nspins)
%
%	Show all (n+1) EPG states in 3x(n+1) plot
%
%	FZ = 3xn column vector of F+,F- and Z.
%	frac = fraction of twist, if animating.
%	scale = axis scale [-scale scale] defaults to 1.
%	Nspins = number of spins, defaults to 24
%
%	Run with no arguments for an example.
%
%	Use arrow3D for prettier plots!

if (nargin < 1) FZ = [.3 0.5; 0.3 0.25; 0.2 0.2]; end;
if (nargin < 2) frac = 0; end;
if (nargin < 3) scale = 1; end;
if (nargin < 4) Nspins = 24; end;


[m,n] = size(FZ); 	% # subplots = size of FZ matrix.	


clf;			% Clear figure
for mm = 1:m
 for nn=1:n
  subplot(m,n,nn+(mm-1)*n)
   Q = 0*FZ; Q(mm,nn)=FZ(mm,nn);		% Just F0 state.

   % -- Setup Titles
   if (mm==1) tt = sprintf('F_{%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==2) tt = sprintf('F_{-%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==3) tt = sprintf('Z_{%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==2) && (nn==1) Q=FZ; tt='All'; end;	% All states in F0* position

   % -- Plot this state
   epg_showstate(Q,frac,scale,Nspins);
   title(tt);
 
   % -- Show F states from above. 
   if  (mm==1)  view(0,90); end;
   if ((mm==2) && (nn>1)) view(0,90); end;
 end;
end;

setprops;
 




