function epg_showstate(FZ,frac,scale,Nspins,voxelvar,simpleaxes)
%function epg_showstate(FZ,frac,scale,Nspins,voxelvar,simpleaxes)
%
%	Show an EPG representation in a single 3D plot.
%	(Zero-out other states if just wanting to show 1 basis function.
%
%	FZ = 3x1 column vector of F+,F- and Z.
%	frac = fraction of twist, if animating.
%	scale = axis scaling [-scale scale] defaults to 1.
%	Nspins = #spins to show, default = 23.
%	voxelvar = origin of spin vectors (0=center, 1=along mz, 2=along mx)
%
%	Get arrow3D to make these look nicer!
%	See epg_show and epg_showorder
%
if (nargin < 4 || length(Nspins)<1) Nspins = 23; end;
myc = mycolors(Nspins);

if (nargin < 1 || length(FZ)<1) FZ = [.75;.25;-.433i]; end;
if (nargin < 2 || length(frac)<1) frac = 0; end;
if (nargin < 3 || length(scale)<1) scale = 1; end;
if (nargin < 4 || length(Nspins)<1) Nspins = 23; end;
if (nargin < 5 || length(voxelvar)<1) voxelvar=0; end;
if (nargin < 6 || length(simpleaxes)<1) simpleaxes=0; end;

% -- Get vectors to plot

M = epg_FZ2spins(FZ,Nspins,frac)*Nspins;

% -- Figure out spin origins (twist vs all at 0 etc)
if (voxelvar==0) spinorig = 0*M; 	% -- All at origin (most common)
else
  spinrange = ([1:Nspins]-.5)/Nspins-.5; 	% voxel dimension	
  spinrange = 2*spinrange;			% fill plots
  if (voxelvar==1) 
    spinorig = [zeros(2,Nspins); spinrange];	% Variation along mz
  else 
    spinorig = [spinrange; zeros(2,Nspins)];	% Variation along mx (for Z_n)
  end;
end;

% -- Plot vectors.
hold off;

for k=1:Nspins
  if (exist('arrow3D'))
    %arrow3D(spinorig(:,k),M(:,k),myc(k,:),0.8,0.03*scale);  % Pre 2021 arrow3D
    if (max(abs(M(:,k))) > 0.01)			
	% -- new arrow3D doesn't like zero-length vectors
	arrow3D(spinorig(:,k),M(:,k),myc(k,:),0.8);
    end;
  else
    Mplot = M+spinorig;		% Add origin, s
    h = plot3([spinorig(1,k) Mplot(1,k)],[spinorig(2,k) Mplot(2,k)],[spinorig(3,k) Mplot(3,k)]);
    set(h,'LineWidth',3);
    set(h,'Color',myc(k,:));
    view(37,30);
    hold on;
    grid on;
  end;
end;
axis(scale*[-1 1 -1 1 -1 1]);

if (voxelvar == 0)
  xlab='M_x';
  ylab='M_y';
  zlab='M_z';
end;
if (voxelvar == 1)
  xlab='M_x';
  ylab='M_y';
  zlab='Voxel Dir (z) x2';
end;
if (voxelvar ==2 )
  xlab='Voxel Dir (z)';
  ylab=' ';
  zlab='M_z';
end;

% Apply labels
xlabel(xlab);
ylabel(ylab);
zlabel(zlab);
   
if (simpleaxes==1)
     hold on;
     axis off;
     plot3([-1 1],[0 0],[0 0]); % x axis
     plot3([0 0],[-1 1 ],[0 0]);        % y axis
     plot3([0 0],[0 0],[-1 1]); % z axis
     text(scale,0,0,xlab); text(0,scale,0,ylab);                      % Mx, My labels
     text(0,0,scale,zlab);
     hold off;
end;

axis square;
lighting phong; camlight right;




function c = mycolors(n)

cc = [
         0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0];

sz = size(cc);
cind = round(([1:n]/(n+1))*sz(1));
c = cc(cind,:);

