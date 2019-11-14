

function epg_show(FZ,frac,scale,Nspins,showtwists)
%function epg_show(FZ,frac,scale,Nspins,showtwists)
%
%	Show all (n+1) EPG states in 3x(n+1) plot
%
%	FZ = 3xn column vector of F+,F- and Z.
%	frac = fraction of twist, if animating.
%	scale = axis scale [-scale scale] defaults to 1.
%	Nspins = number of spins, defaults to 24
%       showtwists = 0 (all spins start at 0) or 1 (variation to show twists)
%
%	Run with no arguments for an example.
%
%	Use arrow3D for prettier plots!

if (nargin < 1 || length(FZ)<1) FZ = [.3 0.5; 0.3 0.25; 0.2 0.2]; end;
if (nargin < 2 || length(frac)<1) frac = 0; end;
if (nargin < 3 || length(scale)<1) scale = 1; end;
if (nargin < 4 || length(Nspins)<1) Nspins = 23; end;
if (nargin < 5 || length(showtwists)<1) showtwists=0; end;


[m,n] = size(FZ); 	% # subplots = size of FZ matrix.	



clf;				% -- Clear figure
for mm = 1:m			% -- Rows (F+,F-,Z)
 for nn=1:n			% -- Order (n)
  subplot(m,n,nn+(mm-1)*n)
   Q = 0*FZ; Q(mm,nn)=FZ(mm,nn);		% Just 1 basis at a time

   % -- Setup Titles
   if (mm==1) tt = sprintf('F_{%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==2) tt = sprintf('F_{-%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==3) tt = sprintf('Z_{%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==2) && (nn==1) Q=FZ; tt='3D View'; end;	% All states in F0* position

   % -- Plot this state
   if (showtwists==0 || (mm==2 && nn==1))	% Don't twist 3D plot
     epg_showstate(Q,frac,scale,Nspins);
   else
     voxvar = 1;		  % -- Variation along mz for F states
     if (mm==3) voxvar = 2; end;  % -- Variation along mx for Z states (easier) 
     epg_showstate(Q,frac,scale,Nspins,voxvar);
   end;
 
   title(tt);
 
   % -- Show F states from above. 
   if  (showtwists==0 && mm==1)  view(0,90); end;
   if ((showtwists==0 && mm==2) && (nn>1)) view(0,90); end;
 end;
end;

setprops;

% -- Save frame if animating.
if (exist('epg_saveframe'))
  epg_saveframe;  % Save frame, if globals 'framenum' and 'filestem' exist
end;

 




