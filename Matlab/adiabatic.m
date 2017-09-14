%
% function [b1,freq,phase] = adiabatic(peakb1,bw,beta,T,Ts,blochsim)
%
%	Function designs a Silver-Hoult adiabatic pulse with
%	peakb1 = max B1 in Gauss.
%	bw = bandwidth in Hz
%	beta = S-H beta (in Hz)
%	T = duration in seconds 
%	Ts = sample period in seconds
%	blochsim = 1 to simulate pulse.
%
%	Run with no arguments for a sample.
%

function [b1,freq,phase] = adiabatic(peakb1,bw,beta,T,Ts,blochsim)

if (nargin < 1) peakb1=0.2;	end;
if (nargin < 2) bw = 2000;	end;
if (nargin < 3) beta = 1000;	end;
if (nargin < 4) T = 0.010;	end;
if (nargin < 5) Ts = 0.00001;	end;
if (nargin < 6) blochsim=0;	end;

T = 2*round(T/Ts/2)*Ts;
N = T/Ts;


t = [Ts:Ts:T] - round(N/2)*Ts;	% Time from -T/2 to T/2.

b1 = peakb1 * sech(beta*t);
freq = bw/2 * tanh(beta*t);
phase = cumsum(freq)*2*pi*Ts;
phase = phase-phase(round(N/2));	% Zero phase half-way.
phase = mod(phase+pi,2*pi)-pi;		% Limit to [-pi,pi];


if (blochsim)
	t = t-t(1);	% Start t at 0 for plots.
	figure(1);
	tt = sprintf('Adiabatic Silver-Hoult Pulse (Beta = %g Hz)',beta);
	subplot(3,1,1);
	plot(t,b1); xlabel('Time(s)'), ylabel('B1(G)');	
	title(tt);
	subplot(3,1,2);
	plot(t,phase); xlabel('Time(s)'), ylabel('Phase(rad)');	
	subplot(3,1,3);
	plot(t,freq); xlabel('Time(s)'), ylabel('Freq(Hz)');	


		
	if (exist('bloch'))
	  gr = 0*b1;
	  tp = Ts;
	  t1 = .6; t2=.1;
	  df = [-3*bw:bw/20:3*bw];
	  dp = 0;
	  mode = 0;
	  [mx,my,mz] = bloch(b1.*exp(i*phase),gr,tp,t1,t2,df,dp,mode);
	
	  figure(2); 
	  plot(df,mx,'b--',df,my,'g--',df,mz,'r-');
	  xlabel('Freq (Hz)');
	  ylabel('Magnetization');
	  legend('Mx','My','Mz');
	  grid;
 	else
	  disp('bloch.m not found.  No Bloch simulation.');
	end;
end;

