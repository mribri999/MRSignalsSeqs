%function showspins(M,scale,spinorig,myc)
%
%	Show vector plot on one axis, that can then be rotated for 2D or 3D
%	viewing.
%
%	M = 3xN spins to show
%	scale = axis scaling [-scale scale] defaults to 1.
%	spinorig = 3xN origin of spins
%	myc = colors to show spins (will default to something reasonable.)
%
%	Get arrow3D to make these look nicer!
%
function showspins(M,scale,spinorig,myc)

if (nargin < 1) M = [0.8,0,0.2].'; end;
if (nargin < 2) scale = 1.0; end;
if (nargin < 3) spinorig = 0*M; end;

sz = size(M);
Nspins = sz(2);
myc = mycolors(Nspins);	% -- Get nice colors, if not passed

% -- Plot vectors.

hold off;
for k=1:Nspins
  if (exist('arrow3D'))
    %arrow3D(spinorig(:,k),M(:,k),myc(k,:),0.8,0.03*scale);
    arrow3D(spinorig(:,k),M(:,k),myc(k,:),0.8);
  else
    h = plot3([spinorig(1,k) M(1,k)],[spinorig(2,k) M(2,k)],[spinorig(3,k) M(3,k)]);
    set(h,'LineWidth',3);
    set(h,'Color',myc(k,:));
    hold on;
    grid on;
  end;
end;
axis(scale*[-1 1 -1 1 -1 1]);

% -- Label Axes, set nice effects
xlabel('M_x');
ylabel('M_y');
zlabel('M_z');
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
cind = ceil(([1:n]/(n+1))*sz(1));
c = cc(cind,:);

