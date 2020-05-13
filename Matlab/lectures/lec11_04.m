%
%	Example of spiral design, sampling, and reconstruction
%
%

% -- Design spiral waveform
FOV = 24;
rmax = 5;	% 5cm^(-1).  1mm resolution.

Fcoeff = FOV;	% Uniform density
N = 120;		% Interleaves
smax = 15000;	% G/cm/s
gmax = 5;	% G/cm
T = .000004;	% 

[k,g,s,time,r,theta] = vds(smax,gmax,T,N,Fcoeff,rmax);


% -- Duplicate Interleaves and calculate density compensation

srot = exp(i*pi*2*[0:N-1]/N);	% Vector to Duplicate and rotate interleaves
sdup = ones(1,N);		% Vector to duplicate	

ksp = k(:)*srot;	% k-space in cm^{-1}	

dcf = vecdcf(g,k);	% Vector DCF appoximation, works well for spiral
dcf = dcf(:)*sdup;	% Duplicate for all interleaves


% -- Generate Sample Data (Array of squares)

[cx,cy] = meshgrid(2*[-2:2],2*[-2:2]);
centers = [cx(:)+i*cy(:); [-10; -8; -6; 6; 8; 10]; i*[-10; -8; -6; 6; 8; 10]]; 
kdata = ksquare(centers,1.9,ksp,T,100);	% k-space data, 100Hz off-resonance
[m,n] = size(kdata);

kdata = kdata + (randn(m,n) + i*randn(m,n))/10;	% Add some noise

% -- Reconstruct Image

kgrid = ksp * (240/256) * (0.5/5);	% Scale to inverse pixels
					% 240cm FOV/ 256 pixels * 
					% 0.5 inv-pixels / 5 cm^(-1)

dat = gridmat(kgrid(:),kdata(:),dcf(:),256);	% Grid data
im = ft(dat);					% FFT to image

dispim(im);			% Display image
title('Reconstructed Image');	% Image

disp('Note the roll-off of the image');
disp('Try changing the noise amplitude');




