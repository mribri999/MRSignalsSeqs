
% Lecture 7, Example 06
% 
% SENSE reconstruction with 4 coils and various R

% -- Start by making a nice image that will show aliasing etc

clear; close all;
N=256;
Nc=4;	% Hard-coded anyway, but cleaner.
R=3;
psi = eye(Nc);	% Equal, uncorrelated noise
cscale = 30000;	% Roughly normalize coil sensitivity

set(0,'defaultAxesFontSize',8); % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width

im1 = 5*circ(N/2,.33*N/2);
im2 = 10*diamond(N/2,.8*N/2);
im3 = 12*diamond(N/2,.5*N/2);
im4 = 10*circ(N/2,.25*N/2);
im = cat(1,cat(2,im1,im2),cat(2,im3,im4));

%dispim(im); title('Truth Image');

% -- Coil Sensitivities
coils=zeros(N,N,4);
coils(:,:,1) = gaussian2d(N,[N/2,N/4],N/4);
coils(:,:,2) = gaussian2d(N,[N/4,N/2],N/4);
coils(:,:,3) = gaussian2d(N,[3*N/4,N/2],N/4);
coils(:,:,4) = gaussian2d(N,[N/2,3*N/4],N/4);
coils = cscale * coils;		% Roughly normalize
figure(1);
% -- Channel images
cims = 0*coils;
for k=1:4 
  cims(:,:,k) = coils(:,:,k).*im; 
  subplot(2,4,k); dispim(coils(:,:,k)); tt = sprintf('Coil Sensitivity %d',k); title(tt); 
  axis equal; axis off;
  subplot(2,4,k+4); dispim(cims(:,:,k)); tt = sprintf('Sig %d',k); title(tt); 
  axis equal; axis off;
end;

% -- Make sure Ny is divisible by R, for simplicity
Ny = floor(N/R)*R;
im = im(1:Ny,:);
coils = coils(1:Ny,:,:);
cims = cims(1:Ny,:,:);

% -- Make Aliased images and sensitivities.
cims = reshape(cims, Ny/R,R,N,Nc);
coils= reshape(coils,Ny/R,R,N,Nc);

figure(2);
coils = permute(coils,[1,3,4,2]);	% reorder sensitivities
for r=1:R
  for c=1:Nc
    subplot(R,Nc,Nc*(r-1)+c);
    dispim(coils(:,:,c,r)); tt=sprintf('Sens Ch %d, p=%d',c,r); title(tt);
    axis equal; axis off;
  end;
end;

% -- Calculate the weights
disp('Calculating weights');
weights = zeros(Ny/R,N,R,Nc);
gfact = zeros(Ny/R,N,R);
noise = zeros(Ny/R,N,R);
for y=1:Ny/R
  for x=1:N
    [w,g,n] = senseweights(squeeze(coils(y,x,:,:)));
    weights(y,x,:,:)=w;

    if (1==0) && (x==128) && (y==60) 
 	sz = size(squeeze(coils(y,x)));
     	disp(sz);
	disp(w); 
    end;
    gfact(y,x,:)=g;
    noise(y,x,:)=sqrt(diag(n));
  end;
end;

% -- Plot the weights
figure(3);
for r=1:R
  for c=1:Nc
    subplot(R,Nc,Nc*(r-1)+c);
    dispim(weights(:,:,r,c)); tt=sprintf('Weight Ch %d, p=%d',c,r); title(tt);
    axis equal; axis off;
  end;
end;

figure(4);
gfact = permute(gfact,[1,3,2]); gfact = reshape(gfact,Ny,N);
noise = permute(noise,[1,3,2]); noise = reshape(noise,Ny,N);
subplot(1,2,1); dispim(gfact); axis off; title('Calculated g-factor');
subplot(1,2,2); dispim(noise); axis off; title('Calculated Noise');
 

% ========= ACTUALLY DO A RECONSTRUCTION!! =======

Nruns = 100;
cims =  permute(cims,[1,3,4,2]);	% reorder coil images
cims =  sum(cims,4);			% ALIAS IMAGES!!

allsenseims = zeros(Ny,N,Nruns);		% Store all!
% -- Generate noise!

tt=sprintf('Doing %d reconstructions',Nruns); disp(tt);
for n=1:Nruns
  chnoise = mvnrnd(zeros(1,Nc),psi,Ny/R*N)+i*mvnrnd(zeros(1,Nc),psi,Ny/R*N);
  chnoise = reshape(chnoise,Ny/R,N,Nc);

  cnims = cims + chnoise;

  % -- Show noisy, aliased channel images
  if (n==1)	% plot only once!
    figure(5);
    for c=1:Nc
      subplot(3,Nc,c); dispim(cims(:,:,c));
      tt=sprintf('Signal, Ch %d',c); title(tt); axis equal; axis off;
      subplot(3,Nc,c+Nc); dispim(chnoise(:,:,c));
      tt=sprintf('Noise, Ch %d',c); title(tt); axis equal; axis off;
      subplot(3,Nc,c+2*Nc); dispim(cnims(:,:,c));
      tt=sprintf('Sig+Noise, Ch %d',c); title(tt); axis equal; axis off;
    end;
  end;

  % -- Reconstruct images!

  senseims = zeros(Ny/R,R,N);	% Set to easily reorder
  
  for x = 1:N
    for y = 1:Ny/R
      csig = squeeze(cnims(y,x,:));
      w = squeeze(weights(y,x,:,:));
      senseims(y,:,x) = w*csig;
   
    end;
  end;
  
  senseims = reshape(senseims,Ny,N);	% Reshape to full image
  allsenseims(:,:,n)=senseims;		% Store result
end;


figure(6);
subplot(2,4,1);
dispim(senseims,0,30); title('Single Image');
subplot(2,4,2);
dispim(mean(allsenseims,3),0,30); title('Mean Image');
subplot(2,4,3);
dispim(std(real(allsenseims),[],3),0,15); title('Std.dev Image');
subplot(2,4,4);
dispim(noise,0,15); title('Calculated. Noise');

x=1:N;
subplot(2,4,5);
plot(x,squeeze(real(allsenseims(floor(Ny/4),:,1))));
lplot('x','signal','1 Image Line, Ny/4');

subplot(2,4,6);
plot(x,squeeze(abs(mean(allsenseims(floor(Ny/4),:,:),3))));
lplot('x','signal','Mean, Line, Ny/4');

subplot(2,4,7);
plot(x,squeeze(std(real(allsenseims(floor(Ny/4),:,:)),[],3)));
lplot('x','signal','Measured Noise, Ny/4');

subplot(2,4,8);
plot(x,squeeze(real(noise(floor(Ny/4),:,:))));
lplot('x','signal','Calculated Noise, Ny/4');
  










