% Lecture 11, Example 03
%
% Radial PSF, with golden-angle sampling
%
% -- Sample in inverse pixels (k from -0.5 to 0.5)
N=200;				% samples along readout
Ngrid=256;			% grid (and image) points
kproj = ([1:N]-0.5-N/2)/N;	% k from -.5 to 0.5
kproj = kproj*N/Ngrid;		% scale so PSF fills FOV
kdcf = abs(kproj);		% Density compensation (1/kr)

ksamps = ones(size(kproj));	% Sample all 1s for PSF
	
Ns = 300;			% (need ~300) for full sampling
ga = 2*pi/(sqrt(5)+1);		% golden angle, about 111 deg.
kangs = [0:Ns-1]*ga;		% sample with golden angle increment.

kall = kproj(:)*exp(i*kangs);	% full k-space locations
kdat = ones(size(kall));	% sample all ones
dcfall=kdcf(:)*ones(1,Ns);	% DCFs for gridding

% -- Plot a few combinations

Twin=[25 50 100 200 300];	% temporal window (#spokes)
Nc = length(Twin);		

set(0,'defaultAxesFontSize',14);        % Default font sizes
for t=1:length(Twin)
  tw=[1:Twin(t)];
  gdat = gridmat(kall(:,tw),kdat(:,tw),dcfall(:,tw),Ngrid);	% grid subset
  psf = ft(gdat);		% PSF is simple reconstruction
  subplot(Nc,2,2*t-1);
  dispim(gdat);	axis off; 
  tt=sprintf('%d Spokes, k-space',Twin(t)); title(tt);
  subplot(Nc,2,2*t); 
  dispim(log(1+abs(psf))); axis off; 			%-- PSF, log scale
  tt=sprintf('%d Spokes, PSF',Twin(t)); title(tt);	  
end;

	


