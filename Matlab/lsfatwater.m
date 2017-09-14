% 	Least-Squares Fat/Water Separation
%
%	dat is organized as kx,ky,TE

dat = fftshift(fft(fftshift(dat,1),[],1),1);	% FFT in x
dat = fftshift(fft(fftshift(dat,2),[],2),2);	% FFT in y
x = size(dat,1);
y = size(dat,2);

Ne = 2;			% Number of echoes to use (1 to Ne)
dTE = 0.5;		% ms (known from the acquisition)
dfcs = .440;		% kHz, fat/water chemical shift frequency

Sw = ones(1,Ne);			% Signal phases for water
Sf = exp(2*pi*i*dfcs*[1:Ne]*dTE);	% Signal phases for fat.

df = [-400:10:400]/1000;		% kHz, range of B0 freq to search
Ndf=length(df);

alldat = zeros(x,y,length(df));		% Allocate, fat and water
allfat = zeros(x,y,length(df));		% Allocate, fat
allwat = zeros(x,y,length(df));		% Allocate, water
allres = zeros(x,y,length(df));		% Allocate, for residuals

dat = reshape(dat(:,:,1:Ne),x*y,Ne);	% Make x*y x Ne
dat = dat.';				% Transpose (conjugate)

for k=1:Ndf
  tt=sprintf('Freq %d of %d',k,Ndf); disp(tt); 	% Display progress

  dfph = exp(2*pi*i*df(k)*[1:Ne]*dTE);	% Phases due to B0 field at TEs
  Smat = diag(dfph)*[Sw(:) Sf(:)];	% Signal matrix:  Si = Smat*[F;W]
  Imat = pinv(Smat);			% pseudoinverse:  [F;W] ~ Imat*Si
  Rmat = Smat*Imat-eye(Ne);		% residual computation:  Si-Si_est

  WFest = Imat * dat;			% Estimate of water & fat
  Res = sum((abs(Rmat * dat)).^2);	% Sum of squares of Act-Fitted Sigs
					% at each TE

  allres(:,:,k) = reshape(Res,x,y);		% Save residual
  allwat(:,:,k) = reshape([1 0]*WFest,x,y);	% Save water (complex)
  allfat(:,:,k) = reshape([0 1]*WFest,x,y);	% Save fat (complex)
  alldat(:,:,k) = reshape([1 i]*abs(WFest),x,y);% Save water & fat mags as I/Q

end;

% -- Try to find the minimum residual over all B0 fields
%	at each pixel

allres = reshape(allres,x*y,Ndf);
allres = allres.';
[yy,m] = min(allres);
allres = allres.';
allres = reshape(allres,x,y,Ndf);
resim = reshape(yy,x,y);			% Keep minimum values of resid.

% -- Extract the appropriate field map values from minimum
fmap = reshape(df(m),x,y)*1000;	% Convert from kHz to Hz.
F=0*fmap;			% Allocate space for fat image
W=F;				% Allocate water image

% -- Extract the water and fat images at each pixel
allwat = reshape(allwat,x*y,Ndf);		% Reshape for extraction
allfat = reshape(allfat,x*y,Ndf);

for k=1:length(m)
  W(k)=allwat(k,m(k));
  F(k)=allfat(k,m(k));
end;

allwat = reshape(allwat,x,y,Ndf);		% Reshape back
allfat = reshape(allfat,x,y,Ndf);

maxf = max(abs(F(:)));
F=F/maxf;					% Scale to max (assume F>W)
W=W/maxf;					% Scale same as fat
fmap = (fmap+max(df)*500)/(max(df)*1000);	% Scale 0-1 (approx)
resim = resim/max(abs(resim(:)));			% Scale 0-1

ims = cat(2,cat(1,W,F),cat(1,fmap,resim));

dispim(ims);				% Crude display, may have to adjust
					% to see Field map, which is in Hz
					% and both positive and negative.
					
  

  


