%function epg_show(FZ,frac,scale,Nspins,showtwists,shiftfs,simpleaxes)
%
%	Show all (n+1) EPG states in 3x(n+1) plot
%
%	FZ = 3xn column vector of F+,F- and Z.
%	frac = fraction of twist, if animating.
%	scale = axis scale [-scale scale] defaults to 1.
%	Nspins = number of spins, defaults to 24
%       showtwists = 0 (all spins start at 0) or 1 (variation to show twists)
%	simpleaxis = 0 (full axes) or 1 for simple, no grid axes%
%	Run with no arguments for an example.
%
%	Use arrow3D for prettier plots!
function epg_show(FZ,frac,scale,Nspins,showtwists,shiftfs,simpleaxes)

if (nargin < 1 || length(FZ)<1) FZ = [.3 0.5; 0.3 0.25; 0.2 0.2]; end;
if (nargin < 2 || length(frac)<1) frac = 0; end;
if (nargin < 3 || length(scale)<1) scale = 1; end;
if (nargin < 4 || length(Nspins)<1) Nspins = 23; end;
if (nargin < 5 || length(showtwists)<1) showtwists=0; end;
if (nargin < 6 || length(shiftfs)<1) shiftfs=1; end;
if (nargin < 7 || length(simpleaxes)<1) simpleaxes=1; end;

[m,n] = size(FZ); 		% # subplots = size of FZ matrix.	
if (m<3) FZ(3,1)=0; end; 	% -- Add 3rd row if not given
if (length(FZ(:)) < 4) shiftfs=0; end;	% No shifting F states if only 0 order




clf;				% -- Clear figure
for mm = 1:m			% -- Rows (F+,F-,Z)
 for nn=1:n			% -- Order (n)
  h=subplotaxis(m,n,nn+(mm-1)*n);

 
  % -- Do plotting for given basis function
   Q = 0*FZ; Q(mm,nn)=FZ(mm,nn);		% Just 1 basis at a time

   % -- Setup Titles
   if (mm==1) tt = sprintf('F{^+_%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
		imag(sum(Q(:)))); end;
   if (mm==2) tt = sprintf('F{^-_%d} = %0.2f + %0.2f',nn-1,real(sum(Q(:))), ...
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

%   Temporary!
    %if (mm==2 && nn==1) axis(4*[-1 1 -1 1 -1 1]); end;
 
   % -- Put simple axes lines instad of grid/white area 
   if (simpleaxes==1) 
     hold on;
     axis off;
     plot3([-1 1],[0 0],[0 0]);	% x axis
     plot3([0 0],[-1 1 ],[0 0]);	% y axis
     plot3([0 0],[0 0],[-1 1]);	% z axis
     text(scale,0,0,'M_x'); text(0,scale,0,'M_y'); 		% Mx, My labels
     if (mm>2 || (mm==2 && nn==1)) text(0,0,scale,'M_z');  end;	% Mz Labels
     hold off;
   end;
 end;
end;

setprops;



% -- Remaining code to shift states if animating - probably ignore!
 
if (shiftfs==1)
 for mm = 1:m			% -- Rows (F+,F-,Z)
  for nn=1:n			% -- Order (n)
   h=subplotaxis(m,n,nn+(mm-1)*n);
   sposx(mm,nn)=h.Position(1);
   sposy(mm,nn)=h.Position(2);
  end;
 end;
 dx = sposx(1,2)-sposx(1,1);	% -- Position between adjacent F's on plot
 dy = sposy(2,1)-sposy(1,1);	% -- Position between rows on plot.

 sposx(1,:) = sposx(1,:)+frac*dx;
 sposx(2,:) = sposx(2,:)-frac*dx;
 sposy(2,2) = sposy(2,2)-frac*dy;	% F(-1) gets shifted left AND up.

 for nn=1:n
   mm=1;
   nnp = n-nn+1;			% Reverse order for 1st row
   h=subplotaxis(m,n,nnp+(mm-1)*n);
   hpos = h.Position;
   set(h,'Position',[sposx(mm,nnp) sposy(mm,nnp) hpos(3:4)]);
 end;

 for nn=1:n
   mm=2;
   nnp = nn;				% Increasing order for 2nd row
   h=subplotaxis(m,n,nnp+(mm-1)*n);
   hpos = h.Position;
   if (nn>1) set(h,'Position',[sposx(mm,nnp) sposy(mm,nnp) hpos(3:4)]); end;
 end;
end;



% -- Save frame if animating.
if (exist('epg_saveframe'))
  epg_saveframe;  % Save frame, if globals 'framenum' and 'filestem' exist
end;


function h=subplotaxis(m,n,q)
% Function calles subplot or subaxis (if it exists), returning the handle.
% subaxis is a more efficient (dense) arrangment of subplots - see Mathworks
% for download.
if (exist('subaxis'))
  h=subaxis(m,n,q);
else
  h=subplot(m,n,q);
end;


