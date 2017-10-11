
%
%	Example 1D SENSE reconstruction with noise.
%

% -- 1. SETUP PARAMETERS	
R=2;				% Reduction Factor

% Coil noise covariance elements 
n12=.3; n13=.13*i;
n14=.1;
n23=.32;
n24=.15;
n34=.28;
ndiag = diag([1 .9 1.2 1]);
nscale = 0.1;			% Noise scaling

cwid = 25;			% Gaussian coil sensitivity width


% -- 2. IMAGE SHAPE 
x = [-64:63];
m =[zeros(1,10) i*[1:20] i*[19:-1:0] zeros(1,10) ...
	sqrt(900-[-30:30].^2) zeros(1,7)];	% 2 shapes

askip = floor(length(x)/R);
x=x(1:askip*R);
m=m(1:askip*R);

% -- 3. COIL SENSITVITIES (Gaussian)
c1 = pdf('norm',x,-45,cwid);	% Coil 1
c2 = pdf('norm',x,-15,cwid);	% Coil 2
c3 = pdf('norm',x,15,cwid);	% Coil 3
c4 = pdf('norm',x,45,cwid);	% Coil 4

c4 = 0*c4;

if (1==0)	% Sanity check!  Uniform coil
  c1 = ones(size(x))/40;
  c2 = 0*c1;
  c3 = 0*c2;
  c4 = 0*c3;
end;

c = 35*[c1(:) c2(:) c3(:) c4(:)]';	% Coil sensitivity matrix
c = diag([1 1.1 0.9 1])*c;		% Different Scaling per coil


% -- 4. CHANNEL IMAGES and K-SPACE
cim = [1;1;1;1]*m .*c;				% Rows of Image * Coil
ksp = ifftshift(ifft(ifftshift(cim,2),[],2),2);	% k-space


% -- 5. DISPLAY STUFF

figure(1);
plotc(x,m);	lplot('x','signal','Image');

figure(2);
plot(x,c);	lplot('x','Coils','Coil Sensitivities');
legend('Coil 1','Coil 2','Coil 3','Coil 4');

figure(3);
plot(x,abs(cim));	lplot('x','Coils','Coil Images');
title('Noise Covariance Matrix \Psi');



% -- 6. GENERATE NOISE (Covariance, then Multivariate Gaussian)

psi = [0 n12 n13 n14; 0 0 n23 n24; 0 0 0 n34; 0 0 0 0];
psi = (psi+psi')+ndiag;
dispim(psi);
n = nscale*(mvnrnd([0 0 0 0],psi,length(x)).' + i*mvnrnd([0 0 0 0],psi,length(x)).');
ksp = ksp + n;


% -- 7. SENSE Reconstruction (SETUP)
ims = 0*x;			% Allocate SENSE image
gfact=ims;			% Allocate g-factor map
gmc = gfact;			% Allocate
gcond = gmc;			% Allocate
inoise = gmc;			% Allocate, for image noise.

askip = floor(length(x)/R);		% Amount to skip for aliased pixel(s).

ipsi = inv(psi);		% Just do inverse-Psi once!
kr = 0*ksp;			% Allocate reduced-k-space matrix
kr(:,1:R:end) = ksp(:,1:R:end);	% Subsampled k-space

ima = R*fftshift(fft(fftshift(kr,2),[],2),2);	% FT to aliased image data.
Cpsave = zeros(R,4,length(x));			% Save Cp values

% -- 8. SENSE (Pixel-wise Recon)

for xx = 1:askip	% Loop across pixels

  C = c(:,xx:askip:end);		% Make Nc x R sensitivity matrix.
  CpC = C'*ipsi*C;			% First step
  iCpC = inv(CpC);			% Second step
  Cp = iCpC*C'*ipsi;			% P-inv Reconstruction matrix.
  Cpsave(:,:,xx) = Cp;			% Save for noise sim.

  dCpC = diag(CpC); 			% Diagonal elements
  diCpC = diag(iCpC);
  
  gfact(xx:askip:end) = sqrt(dCpC.' .* diCpC.');	% gfactor
  gmc(xx:askip:end) = sqrt(dCpC.');	% gfactor - due to coils
  gcond(xx:askip:end) = sqrt(diCpC.');	% gfactor - due to conditioning

  ims(xx:askip:end) = (Cp*ima(:,xx)).';	   % Reconstruct / unalias
  inoise(xx:askip:end) = sqrt(diCpC) * nscale * sqrt(length(x));	% Store calculated noise.


  %tt=sprintf('Pixel Group %d of %d',xx,askip); disp(tt);
end;


% -- 9.  NOISE SIMULATION

Na = 5000;
noiseim = zeros(Na,length(x));

for k = 1:Na		% Noise simulation (per channel)
  %tt=sprintf('Noise run %d of %d',k,Na); disp(tt);
  npts = ceil(length(x)/R);
  cnoise = nscale*sqrt(length(x)) * ...
	(mvnrnd([0 0 0 0],psi,npts).' + i*mvnrnd([0 0 0 0],psi,npts).');
  for xx = 1:askip	% Loop across pixels
    noiseim(k,xx:askip:end) = squeeze(Cpsave(:,:,xx))*cnoise(:,xx);
  end;
end; 
nmean = mean(noiseim);
nstd = std(real(noiseim));



figure(6);
subplot(3,1,1);
plot(x,abs(ima.'));	lplot('x','Coils','Aliased Images (Magnitude)'); 
legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
subplot(3,1,2);
plot(x,real(ima.'));	lplot('x','Coils','Aliased Images (Real)');
legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
subplot(3,1,3);
plot(x,imag(ima.'));	lplot('x','Coils','Aliased Images (Imaginary)'); 
legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
setprops;

figure(7);
plot(x,real(ims),'b--',x,imag(ims),'r--',x,abs(ims),'k-');	
legend('Real','Imag','Magnitude');
lplot('x','Signal','SENSE Image');
setprops;

figure(8);
plot(x,real(gfact),'k-',x,real(gmc),'r--',x,real(gcond),'b--');
legend('g-factor','g_{multicoil}','g_{condition}');
lplot('x','g-factor','g-factor across Image',[min(x) max(x) 0 3]);
setprops;

 	 
figure(9);
plot(x,real(nmean),'b-',x,real(nstd),'r--');
legend('Mean','Std dev');
lplot('x','Noise mean / std-dev','Pseudomultiple Noise Statistics');
setprops;

figure(10);
plot(x,real(inoise));
a = axis; a(3)=-.5; 
lplot('x','Noise','Calculated Noise sqrt(C Psi C)',a);
title('Calculated Noise C^H \Psi C');
setprops;






