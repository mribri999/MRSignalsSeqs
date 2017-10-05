function epg_showorder(FZ,n,frac)
%function epg_showorder(FZ,n,frac)
%
%	Show a single order (n)EPG state in 2x2 plot
%
%	FZ = 3x1 column vector of F+,F- and Z.
%	n = order of state (# twists)
%	frac = fraction of twist, if animating.
%
Nspins = 24;
myc = mycolors(Nspins);

if (nargin < 1) FZ = [.75;.25;-.433i]; end;
if (nargin < 2) n = 1; end;
if (nargin < 3) frac = 0; end;


Q = zeros(3,n+1);
Q(:,n+1)=FZ;

clf;
subplot(2,2,1);
epg_showstate(Q,frac);


subplot(2,2,2);
Q1 = Q; Q1(1:2,:)=0;
epg_showstate(Q1,frac);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
tt = sprintf('Z_{%d} = %0.2f + i%0.2f',n,real(FZ(3)),imag(FZ(3)));
title(tt);


subplot(2,2,3);
Q1 = Q; Q1(2:3,:)=0;
epg_showstate(Q1,frac);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
view(0,90);
tt = sprintf('F_{%d} = %0.2f + i%0.2f',n,real(FZ(1)),imag(FZ(1)));
title(tt);


subplot(2,2,4);
Q1 = Q; Q1(3,:)=0; Q1(1,:)=0;
epg_showstate(Q1,frac);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
view(0,90);
tt = sprintf('F_{-%d} = %0.2f + i%0.2f',n,real(FZ(2)),imag(FZ(2)));
title(tt);
setprops;



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

