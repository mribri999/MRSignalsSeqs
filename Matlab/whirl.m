%
%	function [g] = whirl(N,res,fov)
%
%	Function designs a WHIRL trajectory, as presented by J. Pipe.
%	WHIRL should be slightly faster than spiral and TWIRL, as the
%	maximum gradient amplitude is reached more quickly.  
%
%	INPUT:
%		N	=	# interleaves.
%		res	= 	resolution in cm (NOT mm).
%		fov	= 	field-of-view in cm.
%		dT	= 	Sampling time (default=.000004s)
%
%	OUTPUT:
%		g	= 	gradient waveform.
%
%
%
%	B. Hargreaves,	Oct 2002.
%

%	Note 1:
%		This is sensitive to dT, which bothers me.  Sampling
%		at higher rates, seems to work, so we use an oversampling
%		factor here of 8.

%	See James Pipe, MRM 42(4):714-720, October 1999.

%
%
% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: whirl.m,v $
%	Revision 1.2  2003/02/11 19:41:55  brian
%	Added a bunch of signal simulation files.
%	Not quite sure what the changes to other
%	files were here.
%	
%	Revision 1.1  2002/10/07 22:20:45  brian
%	WHIRL:  Winding Hybrid Interleaved Lines.
%	
%
% ===========================================================



function [g,g1,g2,g3] = whirl(N,res,fov,Ts)

if (nargin < 4)
	Ts=.000004;
end;
upsamp=16;			% Upsample during design.
gmax = 3.9;			% Max grad amp, G/cm.
S = 14500;			% Max slew rate, G/cm/s
dT = Ts/upsamp;			% Sampling rate, seconds
kmax = .5/res;			% maximum k-space radius.

gamma = 4258;			% gamma, in rads/G.
delta = N/2/pi/fov;		% delta, as defined in Pipe.
r1 = 1/sqrt(5)*delta;		% r1, end-radius of stage 1 as in Pipe.
r2 = 3/sqrt(5)*delta;		% r2, end-radius of stage 2 as in Pipe.


tt = sprintf('Delta = %g cm^(-1)',delta); disp(tt);



% =========================== Stage 1 ============================

Gc = sqrt(2*S*r1/gamma);	% Pipe, eq. 26.
g = [0:S*dT:Gc].';		% Define g during ramp.
k = cumsum(g)*gamma*dT;		% Define k during ramp.
ng = length(g);			% Number of points in ramp.
ng1=ng;

tt=sprintf('End stage 1:  G = %g G/cm and kr= %g cm^(-1)',max(abs(g(:))),max(abs(k)));
disp(tt);

% =========================== Stage 2 ============================

%  --- G, r, kk and phi for stages 2 and 3, with values at end of stage 1.
G = g(end);	% G(t) as defined in Pipe.
r = k(end);	% r(t), k-space radius as defined in Pipe.
kk = k(end);	% kk is just the current value of k(t) in Pipe paper.
phi = 0;	% phi(t) as defined in Pipe.

%  --- Pre-allocate memory to speed up the loop ----

maxng = 10000*upsamp;
g(maxng)=0;
k(maxng)=0;
Grec=g;
Grec(10000)=0;
phirec = 0*Grec;

%  --- Circular bridge between trajectories ----
%
%	G(t) is constant.  phi'(t) is constant.
%
done=0;	
dphi = S/abs(G)*dT;	% Pipe, eq. 15.  We use the G that we got to in
			%		 Stage 1 which should be less than,
			% 		 but approach the Gc in Pipe.

while (r < r2)
	ng = ng+1;			% update # gradient points.
	phi = phi+dphi;			% update phi.
	g(ng) = G*exp(i*phi);		% update g(t).
	kk = kk+g(ng)*dT*gamma;		% update k
	k(ng) = kk;			% record k(t).
	r = abs(kk);			% update k-space radius.

end;
ng2=ng;

tt=sprintf('End stage 2:  G = %g G/cm and kr= %g cm^(-1)',max(abs(g(:))),max(abs(k)));
disp(tt);
	
% =========================== Stage 3 ============================

dG=0;
%maxng = ng;	% STOP, if testing!!

while (r<kmax) & (ng < maxng)

	%	G'(t) and phi'(t) are defined in terms of each other,
	%	so we'll take a guess at r and G in the middle of the
	%	sample interval, extrapolating from the last sample.

	Gapp=G;					% approx. G in center.
	r=abs(k(ng) + gamma*g(ng)*dT/2);	% approx. r in center (using k and G).

	ng = ng+1;				% update # gradient points.

	%	dphi = phi'(t)*dT
	dphi = gamma*Gapp/sqrt(r^2-delta^2)*dT;	% Pipe, eq. 14.
	dG = sqrt(S^2-(dphi/dT*Gapp)^2)*dT;	% Pipe, eq. 27.

	G = G+dG;				% update G
	if (G>=gmax)
		G=gmax;				% constrain G.
		dG=0;				% for approx. G above.
	end;

	phi = phi+dphi;				% update phi.
	g(ng) = G*exp(i*phi);			% update g(t).  Pipe Eq 2
	kk = kk+g(ng)*dT*gamma;			% update k
	k(ng) = kk;				% record k(t).
	r = abs(kk);				% update k-space radius.

	if (ng >= maxng)
		tt=sprintf('Maximum number of points, %d, reached.',maxng);
		disp(tt);
	end;

end;
ng3=ng;


% ======================== Finished! ============================

k = k(1:upsamp:ng);
g = g(1:upsamp:ng);
t = [1:length(g)].'*Ts;

g1=g(1:round(ng1/upsamp));
g2=g(round(ng1/upsamp)+1:round(ng2/upsamp));
g3=g(round(ng2/upsamp)+1:length(g));

%  ------ Plot, for testing. -----------
if (1==1)
    figure(1);
    plotgradinfo(g(:),Ts);


    figure(2);
    for q=1:N
	plot(k*exp(2*pi*i*q/N));
	hold on;
    end;
    hold off;
    axis equal;
    title('K-space trajectory...');
end;


