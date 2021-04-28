% function [s,phasediag,FZ] = epg_gradecho(flipangle,T1,T2,TR,TE,rfdphase,rfdphaseinc,gspoil,dfreq,N)
%
%	EPG simulation of any gradient-echo sequence:
%	  (a) balanced SSFP 
%	  (b) gradient spoiling
%	  (c) reversed gradient spoiling
%	  (d) RF spoiling
%	  (e) Perfect spoiling
%
%	INPUT:
%	  flipangle = flip angle in radians [pi/6 radians]
%	  T1,T2 = relaxation times, ms [1000,200 ms]
%	  TR = repetition time, ms [10ms]
%	  TE = echo time (RF to signal), ms [5ms]
%         rfdphase = amount to increment RF phase per TR, radians [pi radians]
%         rfdphaseinc = amount to increment rfdphase, per TR, radians [0 rad]
%         gspoil = 0 to apply no spoiler, 1,-1 spoiler after,before TE [0]
%	  dfreq = off-resonance frequency (Hz) [0Hz]
%	  N = number of repetitions to steady state [200]
%
%	NOTE SETTINGS FOR SEQUENCES:
%	  SEQUENCE             rfdphase            rfdphaseinc    gspoil
%         --===---------------------------------------------------------------
%	 (a) bSSFP           any (shifts profile)         0          0
%        (b) Grad-spoiled          anything               0          1
%        (c) Rev.Grad-spoiled      anything               0          -1
%        (d) RF Spoiled            anything	    pi/180*117       1	
%	 (e) Perfect spoiling	   anything         0 or any         100
%
%	All states are kept, for demo purposes, though this 
%	is not really necessary.  
%
%       OUTPUT:
%	  s = signal
%	  phasediag = display-able phase diagram
%	  FZ = record of all states (3xNxN)

function [s,phasediag,FZ] = epg_gradecho(flipangle,T1,T2,TR,TE,rfdphase,rfdphaseinc,gspoil,dfreq,N)

if (nargin < 1) flipangle = pi/6; end;
if (nargin < 2) T1 = 1000; end; %  ms
if (nargin < 3) T2 = 200; end; %ms 
if (nargin < 4) TR = 10; end;	% ms
if (nargin < 5) TE = 5; end;	% ms
if (nargin < 6) rfdphase = pi; end;  % radians
if (nargin < 7) rfdphaseinc = 0; end;  % radians
if (nargin < 8) gspoil = 0; end;
if (nargin < 9) dfreq = 0; end;
if (nargin < 10) N=200; end;

% -- Initialization
P = zeros(3,N);			% Allocate all known states, 1 per echo.
P(3,1)=1;			% Initial condition/equilibrium.
FZ = zeros(3,N,N);		% Store all states.

Pstore = zeros(2*N,N);		% Store F,F* states, with 2 gradients per echo
Zstore = zeros(N,N); 		% Store Z states, with 2 gradients per echo

s = zeros(1,N);		% Allocate signal vector to store.
rfphase = 0;		% Register to store.

for n=1:N	% Repeat over N TRs.  *START/END at TE!!*

  P = epg_grelax(P,T1,T2,TR-TE,0,0,0,1);   % Relaxation TE to TR
  P = epg_zrot(P,2*pi*dfreq*(TR-TE)/1000);  % Off-resonance precession TE to TR
  if (gspoil ==1) P = epg_grad(P,1); end;  % Gradient spoiler at end of TR
  if (gspoil==100) P = diag([0 0 1])*P; end;% Zero-out transverse states at TR

  P = epg_rf(P,flipangle,rfphase);	% RF excitation

  P = epg_grelax(P,T1,T2,TE,0,0,0,1);   % Relaxation 0 to TE
  P = epg_zrot(P,2*pi*dfreq*(TE)/1000); % Off-resonance precession 0 to TE.
  if (gspoil < 0) P = epg_grad(P,1); end;  % Gradient spoiler before TE

  % -- Store Signals and update RF phase
  s(1,n) = P(1,1)*exp(-i*rfphase);	% Store signal, demodulated by RF phase
  rfphase = rfphase+rfdphase;		% Increment the phase
  rfdphase = rfdphase + rfdphaseinc;	% Increment the increment...

  Pstore(N:2*N-1,n) = P(2,:).';		% Put in negative states
  Pstore(1:N,n) = flipud(P(1,:).');  	% Put in positive, overwrite center.
  Zstore(:,n) = P(3,:).';
  FZ(:,:,n) = P;			% Store in state matrix format
end;

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width

plotstate = cat(1,Pstore,Zstore);
subplot(1,2,1);
plot([1:N]*TR,abs(s));
lplot('Evolution Time','Signal',' Signal vs TR time');

subplot(1,2,2);
dispim(plotstate,0,0.2);
xlabel('Echo'); ylabel('F(top) and Z(bottom) States');
title('F and Z states vs time');
phasediag = plotstate;	

