%
%	function [k,g,s,time,r,theta] = vds(smax,gmax,T,N,Fcoeff,rmax,z)
%
%
%	VARIABLE DENSITY SPIRAL GENERATION:
%	----------------------------------
%
%	Function generates variable density spiral which traces
%	out the trajectory
%				 
%			k(t) = r(t) exp(i*q(t)), 		[1]
%
%	Where q is the same as theta...
%		r and q are chosen to satisfy:
%
%		1) Maximum gradient amplitudes and slew rates.
%		2) Maximum gradient due to FOV, where FOV can
%		   vary with k-space radius r/rmax, as
%
%			FOV(r) = Sum    Fcoeff(k)*(r/rmax)^(k-1)   [2]
%
%
%	INPUTS:
%	-------
%	smax = maximum slew rate G/cm/s
%	gmax = maximum gradient G/cm (limited by Gmax or FOV)
%	T = sampling period (s) for gradient AND acquisition.
%	N = number of interleaves.
%	Fcoeff = FOV coefficients with respect to r - see above.
%	rmax= value of k-space radius at which to stop (cm^-1).
%		rmax = 1/(2*resolution);
%	z = R/L for gradient coil, to include voltage model.  This
%		tends to speed up gradient design slightly.
%
%
%	OUTPUTS:
%	--------
%	k = k-space trajectory (kx+iky) in cm-1.
%	g = gradient waveform (Gx+iGy) in G/cm.
%	s = derivative of g (Sx+iSy) in G/cm/s.
%	time = time points corresponding to above (s).
%	r = k-space radius vs time (used to design spiral)
%	theta = atan2(ky,kx) = k-space angle vs time.
%
%
%	METHODS:
%	--------
%	Let r1 and r2 be the first derivatives of r in [1].	
%	Let q1 and q2 be the first derivatives of theta in [1].	
%	Also, r0 = r, and q0 = theta - sometimes both are used.
%	F = F(r) defined by Fcoeff.
%
%	Differentiating [1], we can get G = a(r0,r1,q0,q1,F)	
%	and differentiating again, we get S = b(r0,r1,r2,q0,q1,q2,F)
%
%	(functions a() and b() are reasonably easy to obtain.)
%
%	FOV limits put a constraint between r and q:
%
%		dr/dq = N/(2*pi*F)				[3]	
%
%	We can use [3] and the chain rule to give 
%
%		q1 = 2*pi*F/N * r1				[4]
%
%	and
%
%		q2 = 2*pi/N*dF/dr*r1^2 + 2*pi*F/N*r2		[5]
%
%
%
%	Now using [4] and [5], we can substitute for q1 and q2
%	in functions a() and b(), giving
%
%		G = c(r0,r1,F)
%	and 	S = d(r0,r1,r2,F,dF/dr)
%
%
%	Using the fact that the spiral should be either limited
%	by amplitude (Gradient or FOV limit) or slew rate, we can
%	solve 
%		|c(r0,r1,F)| = |Gmax|  				[6]
%
%	analytically for r1, or
%	
%	  	|d(r0,r1,r2,F,dF/dr)| = |Smax|	 		[7]
%
%	analytically for r2.
%
%	[7] is a quadratic equation in r2.  The smaller of the 
%	roots is taken, and the real part of the root is used to
%	avoid possible numeric errors - the roots should be real
%	always.
%
%	The choice of whether or not to use [6] or [7], and the
%	solving for r2 or r1 is done by findq2r2 - in this .m file.
%
%	Once the second derivative of theta(q) or r is obtained,
%	it can be integrated to give q1 and r1, and then integrated
%	again to give q and r.  The gradient waveforms follow from
%	q and r. 	
%
%	Brian Hargreaves -- Sept 2000.
%
%	See Brian's journal, Vol 6, P.24.
%
%
%	See also:  vds2.m,  vdsmex.m,  vds.c
%

% =============== CVS Log Messages ==========================
%	$Log: vds.m,v $
%	Revision 1.5  2004/04/27 18:08:44  brian
%	Revision 1.6  2007-02-23 18:30:19  brian
%		* Included LR-circuit model for increased slewing.
%		* NOTE:  There is a bug in the substititution of
%			q1 with r1, in that the dr/dt term should
%			include dFOV/dr.  This is confusing, but should
%			be a relatively easy fix that stabilizes the
%			algorithm for quickly-changing density functions.
%
%	Revision 1.5  2004/04/27 18:08:44  brian
%	Changed FOV to a polynomial of unlimited length,
%	and hopefully changed all comments accordingly.
%	Also moved sub-functions into vds.m so that
%	no other .m files are needed.
%	
%	Revision 1.4  2003/09/16 02:55:52  brian
%	minor edits
%	
%	Revision 1.3  2002/11/18 05:36:02  brian
%	Rounds lengths to a multiple of 4 to avoid
%	frame size issues later on.
%	
%	Revision 1.2  2002/11/18 05:32:19  brian
%	minor edits
%	
%	Revision 1.1  2002/03/28 01:03:20  bah
%	Added to CVS
%	
%
% ===========================================================

function [k,g,s,time,r,theta] = vds(smax,gmax,T,N,Fcoeff,rmax,z)

if (nargin < 7)
	z=0;
end;

disp('vds.m');
gamma = 4258;

oversamp = 8;		% Keep this even.
To = T/oversamp;	% To is the period with oversampling.



q0 = 0;	
q1 = 0;
theta = zeros(1,10000);
r = zeros(1,10000);
r0 = 0;
r1 = 0;

time = zeros(1,10000);
t = 0;
count = 1;

theta = zeros(1,1000000);
r = zeros(1,1000000);
time = zeros(1,1000000);

while r0 < rmax
	[q2,r2] = findq2r2(smax,gmax,r0,r1,To,T,N,Fcoeff,rmax,z);

	% Integrate for r, r', theta and theta' 	
	q1 = q1 + q2*To;
	q0 = q0 + q1*To;
 	t = t + To;

	r1 = r1 + r2*To;
	r0 = r0 + r1*To;

	% Store.
	count = count+1; 
	theta(count) = q0;
	r(count) = r0;
	time(count) = t;

	if (rem(count,100)==0)
		tt = sprintf('%d points, |k|=%f',count,r0);
		disp(tt);
	end;
end;

r = r(oversamp/2:oversamp:count);
theta = theta(oversamp/2:oversamp:count);
time = time(oversamp/2:oversamp:count);

%	Keep the length a multiple of 4, to save pain...!
%
ltheta = 4*floor(length(theta)/4);
r=r(1:ltheta);
theta=theta(1:ltheta);
time=time(1:ltheta);

%
% 	Plot.
%
%x = alpha*theta .* cos(theta);
%y = alpha*theta .* sin(theta);

%plot(x,y);
%title('k-space trajectory.');


k = r.*exp(i*theta);

g = 1/gamma*([k 0]-[0 k])/T;
g = g(1:length(k));

s = ([g 0]-[0 g])/T;
s = s(1:length(k));


% ========= Plot gradients and slew rates. ==========

tp= time(1:10:end);
gp = g(1:10:end);		% Plotting undersamples
kp = k(1:10:end);
sp= s(1:10:end);


subplot(2,2,1);
plot(real(kp),imag(kp));
lplot('kx (cm^{-1})','ky (cm^{-1})','ky vs kx');
axis('square');

subplot(2,2,2);
plot(tp,real(kp),'c-',tp,imag(kp),'g-'); %,tp,abs(kp),'k-');
lplot('Time (s)','k (cm^{-1})','k-space vs Time');


subplot(2,2,3);
plot(tp,real(gp),'c-',tp,imag(gp),'g-'); %,tp,abs(gp),'k-');
lplot('Time (s)','g (G/cm)','Gradient vs Time');

subplot(2,2,4);
plot(tp,real(sp),'c-',tp,imag(sp),'g-',tp,abs(sp),'k-');
lplot('Time (s)','s (G/cm/s)','Slew-Rate vs Time');

if (exist('setprops')) setprops; end;

return;




%
%  	function [q2,r2] = q2r2(smax,gmax,r,r1,T,Ts,N,F)
%
%	VARIABLE DENSITY SPIRAL DESIGN ITERATION
%	----------------------------------------
%	Calculates the second derivative of r and q (theta),
%	the slew-limited or FOV-limited
%	r(t) and q(t) waveforms such that 
%
%		k(t) = r(t) exp(i*q(t))
%
%	Where the FOV is a function of k-space radius (r)
%
%	FOV = Fcoeff(1) + Fcoeff(2)*r/rmax + Fcoeff(3)*(r/rmax)^2 + ... ;
%
%	F(1) in cm.
%	F(2) in cm^2.
%	F(3) in cm^3.
%	.
%	.
%	.
%
%	The method used is described in vds.m
%	
%	INPUT:
%	-----
%	smax  	= Maximum slew rate in G/cm/s.
%	gmax 	= Maximum gradient amplitdue in G.
%	r	= Current value of r.
%	r1 	= Current value of r', first derivative of r wrt time.
%	T	= Gradient sample rate.
%	Ts	= Data sampling rate.
%	N	= Number of spiral interleaves.
%	F is described above.
%


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: vds.m,v $
%	Revision 1.5  2004/04/27 18:08:44  brian
%	Changed FOV to a polynomial of unlimited length,
%	and hopefully changed all comments accordingly.
%	Also moved sub-functions into vds.m so that
%	no other .m files are needed.
%	
%	Revision 1.2  2003/05/29 23:02:21  brian
%	minor edits
%	
%	Revision 1.1  2002/03/28 01:03:20  bah
%	Added to CVS
%	
%
% ===========================================================


	
function [q2,r2] = findq2r2(smax,gmax,r,r1,T,Ts,N,Fcoeff,rmax,z)

gamma = 4258;			% Hz/G

smax = smax + z*gmax;

F = 0;			% FOV function value for this r.
dFdr = 0;		% dFOV/dr for this value of r.
for rind = 1:length(Fcoeff)
	F = F+Fcoeff(rind)*(r/rmax)^(rind-1);
	if (rind>1)
		dFdr = dFdr + (rind-1)*Fcoeff(rind)*(r/rmax)^(rind-2)/rmax;
	end;
end;

GmaxFOV = 1/gamma /F/Ts;		% FOV limit on G
Gmax = min(GmaxFOV,gmax);	%

maxr1 = sqrt((gamma*Gmax)^2 / (1+(2*pi*F*r/N)^2));  


if (r1 > maxr1)			
			% Grad amplitude limited.  Here we
			% just run r upward as much as we can without
			% going over the max gradient.
	r2 = (maxr1-r1)/T; 
	%tt = sprintf('Grad-limited r=%5.2f, r1=%f',r,r1);
	%disp(tt);

else

	twopiFoN = 2*pi*F/N;
	twopiFoN2 = twopiFoN^2;

	%	A,B,C are coefficents of the equation which equates
	% 	the slew rate calculated from r,r1,r2 with the
	%	maximum gradient slew rate.
	%
	%	A*r2*r2 + B*r2 + C  =  0	
	%
	%	A,B,C are in terms of F,dF/dr,r,r1, N and smax.
	%

	A = 1+twopiFoN2*r*r;
	B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1 + 2*z*r1 + 2*twopiFoN2*r1*r;;
	C1 = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4 - (gamma)^2*smax^2;
	C2 = z*(z*r1^2 + z*twopiFoN2*r1^2 + 2*twopiFoN2*r1^3*r + 2*twopiFoN2/F*dFdr*r1^3*r);
	C = C1+C2;


	%A = 1+twopiFoN2*r*r;
	%B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1;
	%C = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4 - (gamma)^2*smax^2;


	[rts] = qdf(A,B,C);	% qdf = Quadratic Formula Solution.
	r2 = real(rts(1));	% Use bigger root.  The justification
				% for this is not entirely clear, but
				% in practice it seems to work, and 
				% does NOT work with the other root.


	% Calculate resulting slew rate and print an error 
	% message if it is too large.
	
	slew = 1/gamma*(r2-twopiFoN2*r*r1^2 + i*twopiFoN*(2*r1^2 + r*r2 + dFdr/F*r*r1^2));
	%tt = sprintf('Slew-limited r=%5.2d  SR=%f G/cm/s',r,abs(slew));
	%disp(tt);
	sr = abs(slew)/smax;

	if (abs(slew)/smax > 1.01)
		tt = sprintf('Slew violation, slew = %d, smax = %d, sr=%f, r=%f, r1=%f',round(abs(slew)),round(smax),sr,r,r1);
		disp(tt);
	end;

end;


%	Calculate q2 from other pararmeters.

q2 = 2*pi/N*dFdr*r1^2 + 2*pi*F/N*r2;






%	function [r1,r2] = qdf(a,b,c)
%
%	Outputs quadratic roots of ax^2+bx+c = 0.
%


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: vds.m,v $
%	Revision 1.6  2007-02-23 18:30:19  brian
%		* Included LR-circuit model for increased slewing.
%		* NOTE:  There is a bug in the substititution of
%			q1 with r1, in that the dr/dt term should
%			include dFOV/dr.  This is confusing, but should
%			be a relatively easy fix that stabilizes the
%			algorithm for quickly-changing density functions.
%
%	Revision 1.5  2004/04/27 18:08:44  brian
%	Changed FOV to a polynomial of unlimited length,
%	and hopefully changed all comments accordingly.
%	Also moved sub-functions into vds.m so that
%	no other .m files are needed.
%	
%	Revision 1.1  2002/03/28 01:27:46  bah
%	Added to CVS
%	
%
% ===========================================================


function [roots] = qdf(a,b,c)

d = b^2 - 4*a*c;

roots(1) = (-b + sqrt(d))/(2*a);
roots(2) = (-b - sqrt(d))/(2*a);










