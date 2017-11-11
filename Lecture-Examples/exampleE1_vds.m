%% Variable-density sampling (1D)
%
N = 64;
Ni=100;
denscomp=0;			% Density compensate
kx = [-N/2:1/Ni:N/2-1];		% kx locations
a = .1;				% Weighting for spacing
dk = 1+a*abs(kx);		% Spacing betweek k-space points
H = 0*kx;

k=0;
fc = find(kx==0);

while k<max(kx)
  [y,imin] = min(abs(kx-k));
  H(imin)=1;
  if (denscomp==1) H(imin)=dk(imin)/max(abs(dk)); end;
  H(2*fc-imin)=H(imin);	% Negative point
  k=k+dk(imin);
end;

% -- Calculate PSFs

Hv=H;
sH = sum(Hv);

% -- Full k-space, PSF calculation
Hf = ones(1,N);
H=Hf;
Hpad = [0*H 0*H 0*H 0*H 0*H 0*H H 0*H 0*H 0*H 0*H 0*H 0*H];
hf = ifftshift(ifft(ifftshift(Hpad)))*length(Hpad)/length(Hf);
xf = [0:length(hf)-1]/length(hf)*N-N/2;
hf = hf.* exp(i*pi*xf/64);

% -- VDS k-space, PSF calculation
H = Hv;
Hpad = [0*H 0*H 0*H 0*H 0*H 0*H H 0*H 0*H 0*H 0*H 0*H 0*H];
hv = ifftshift(ifft(ifftshift(Hpad)))*length(Hpad)/length(H)*Ni*N/sH;
xv = (([0:length(hv)-1]-0.5)/length(hv)*N-N/2)*Ni;
hv = hv.* exp(-i*pi*xv/64);

% -- Plots
subplot(2,1,1)
stem(kx,Hv(:));
axis([-N/2 N/2-1 -.1 1.1]);
lplot('k-space','Sample Modulation','Sampling Modulation H(k)');

subplot(2,1,2)
setprops
plot(xf,real(hf),'b--',xv,real(hv),'r--');
lplot('position (pixels)','PSF','PSF h(r)');
%axis([-N/2 N/2-1 -.6 1.2]);
axis([-5 5 -.6 1.2]);
setprops
legend('Full','Var-Dens');

