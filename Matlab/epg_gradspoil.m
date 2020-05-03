% function [s,phasediag,P] = epg_gradspoil(flipangle,N,T1,T2,TR)
%
%	EPG Simulation of Gradient-spoiled sequence.  
%
%	flipangle = flip angle (radians)
%	N = number of TRs
%	T1,T2,TR = relaxation times and TR
%
%	All states are kept, for demo purposes, though this 
%	is not really necessary.

function [s,phasediag,P] = epg_gradspoil(flipangle,N,T1,T2,TR)

if (nargin < 1) flipangle = pi/6; end;
if (nargin < 2) N=100; end;
if (nargin < 3) T1 = 1; end; %  sec
if (nargin < 4) T2 = 0.1; end; % sec
if (nargin < 5) TR = 0.01; end;	% sec


P = zeros(3,N);			% Allocate all known states, 1 per echo.
P(3,1)=1;			% Initial condition/equilibrium.

Pstore = zeros(2*N,N);		% Store F,F* states, with 2 gradients per echo
Zstore = zeros(N,N); 		% Store Z states, with 2 gradients per echo

s = zeros(1,N);		% Allocate signal vector to store.

for n=1:N
  P = epg_rf(P,flipangle,pi/2);		% RF excitation
  s(n) = P(1,1);  			% Signal is F0 state.
  P = epg_grelax(P,T1,T2,TR,1,0,1,1);   % Gradient spoiler, relaxation
  Pstore(N:2*N-1,n) = P(2,:).';		% Put in negative states
  Pstore(1:N,n) = flipud(P(1,:).');  	% Put in positive, overwrite center.
  Zstore(:,n) = P(3,:).';
end;

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width

plotstate = cat(1,Pstore,Zstore);
subplot(1,2,1);
plot([1:N]*TR,abs(s));
lplot('Evolution Time','Signal','Gradient-Spoiled Signal vs TR time');

subplot(1,2,2);
dispim(plotstate,0,0.2);
xlabel('Echo'); ylabel('F(top) and Z(bottom) States');
title('F and Z states vs time');
phasediag = plotstate;	

