function epg_showstate(FZ,frac,scale,Nspins)
%function epg_showstate(FZ,frac,scale,Nspins)
%
%	Show an EPG representation state in 3D plot.
%	(Zero-out other states if just wanting to show 1 state)
%
%	FZ = 3x1 column vector of F+,F- and Z.
%	frac = fraction of twist, if animating.
%	scale = axis scaling [-scale scale] defaults to 1.
%	Nspins = #spins to show, default = 24.
%
%	Get arrow3D to make these look nicer!
%	See epg_show and epg_showorder
%
if (nargin < 4) Nspins = 24; end;
myc = mycolors(Nspins);

if (nargin < 1) FZ = [.75;.25;-.433i]; end;
if (nargin < 2) frac = 0; end;
if (nargin < 3) scale = 1; end;

M = epg_FZ2spins(FZ,Nspins,frac)*Nspins;
hold off;
for k=1:Nspins
  if (exist('arrow3D'))
    arrow3D([0 0 0],M(:,k),myc(k,:),0.8,0.03*scale);
  else
    h = plot3([0 M(1,k)],[0 M(2,k)],[0 M(3,k)]);
    set(h,'LineWidth',3);
    set(h,'Color',myc(k,:));
    hold on;
    grid on;
  end;
end;
axis(scale*[-1 1 -1 1 -1 1]);

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
cind = round(([1:n]/(n+1))*sz(1));
c = cc(cind,:);

