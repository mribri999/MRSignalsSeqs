% Lecture 11, Example 02
%
% Radial PSF
%
% -- Sample in inverse pixels (k from -0.5 to 0.5)
N=200;				% samples along readout
Ngrid=256;			% grid (and image) points
kproj = ([1:N]-0.5-N/2)/N;	% k from -.5 to 0.5
kproj = kproj*N/Ngrid;		% scale so PSF fills FOV
kdcf = abs(kproj);		% Density compensation (1/kr)

ksamps = ones(size(kproj));	% Sample all 1s for PSF
	
Ns = 50;			% 150 spokes (need ~300) for full sampling
kangs = pi*[0:Ns-1]/Ns;		% sample over 180 degrees (pi radians)

kall = kproj(:)*exp(i*kangs);	% full k-space locations
kdat = ones(size(kall));	% sample all ones
dcfall=kdcf(:)*ones(1,Ns);	% DCFs for gridding

gdat = gridmat(kall,kdat,dcfall,Ngrid);	% Grid to 256 grid

psf = ft(gdat);		% PSF is simple reconstruction


set(0,'defaultAxesFontSize',14);        % Default font sizes
subplot(2,2,1);				 
plot(kall);					% Nice plot of k-space	
lplot('k_x','k_y','k-space trajectory');	% Label plot

subplot(2,2,2);					
dispim(gdat);	axis off; title('Gridded ones');  % Display gridded k-locs

subplot(2,2,3);
dispim(log(1+abs(psf))); axis off; title('PSF');  % PSF, log scale

subplot(2,2,4);	
dispim(cropim(psf,16,16));	% central PSF (linear scale);
axis off; title('Central PSF');

	


