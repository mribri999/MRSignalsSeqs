function epg_showorder(FZ,n,frac,Nspins,showtwist)
%function epg_showorder(FZ,n,frac,Nspins,showtwist)
%
%	Show a single order (n)EPG state in 2x2 plot
%
%	FZ = 3x1 column vector of F+,F- and Z. [default .75;.25;-.433i]
%	n = order of state (# twists) [default = 1]
%	frac = fraction of twist, if animating. [default=0]
%	Nspins = number of spins [default=23]
%	showtwist = 1 to show twists.


if (nargin < 1 || length(FZ)<1) FZ = [.75;.25;-.433i]; end;
if (nargin < 2 || length(n)<1) n = 1; end;
if (nargin < 3 || length(frac)<1) frac = 0; end;
if (nargin < 4 || length(Nspins)<1) Nspins = 23; end;
if (nargin < 5 || length(showtwist)<1) showtwist = 0; end;


if (showtwist==1)
  ftwist=1;
  ztwist=2;
else
  ftwist=0;
  ztwist=0;
end;



Q = zeros(3,n+1);
Q(:,n+1)=FZ;

clf;
subplot(2,2,1);
epg_showstate(Q,frac,1,Nspins);
title('3D View')

subplot(2,2,2);
Q1 = Q; Q1(1:2,:)=0;
epg_showstate(Q1,frac,1,Nspins,ztwist);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
tt = sprintf('Z_{%d} = %0.2f + i%0.2f',n,real(FZ(3)),imag(FZ(3)));
title(tt);


subplot(2,2,3);
Q1 = Q; Q1(2:3,:)=0;
epg_showstate(Q1,frac,1,Nspins,ftwist);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
if (showtwist==0) view(0,90); end;
tt = sprintf('F_{%d} = %0.2f + i%0.2f',n,real(FZ(1)),imag(FZ(1)));
title(tt);


subplot(2,2,4);
Q1 = Q; Q1(3,:)=0; Q1(1,:)=0;
epg_showstate(Q1,frac,1,Nspins,ftwist);
axis([-1 1 -1 1 -1 1]);
lighting phong; camlight right;
if (showtwist==0) view(0,90); end;

tt = sprintf('F_{-%d} = %0.2f + i%0.2f',n,real(FZ(2)),imag(FZ(2)));
title(tt);
setprops;

