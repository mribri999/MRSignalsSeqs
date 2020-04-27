
% Lecture 7, Example 04
% 
% Starting with k-space, change NEX, FOV, Resolution to 
% test the impact on SNR
%
% -- Defaults
set(0,'defaultAxesFontSize',8);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width

% -- Generate k-space for a square object with large FOV
N=256;			% Image size (square)
siglevel=10;		% Signal level in square
im = zeros(N,N);
im(N/2-N/8+1:N/2+N/8, N/2-N/8+1:N/2+N/8) = ones(N/4,N/4)*siglevel;
ksp = N*ift(im);	% iFFT, normalized

roi = find(im(:)>0.5);	% Include all points in original image above 0.5
bgroi = [1:8*N];	% First 8 lines, assumed to be in background
 
% -- Add noise, and repeat multiple times
Nreps=500;
kall = zeros(N,N,Nreps);
for n = 1:Nreps kall(:,:,n)=ksp; end;
kall = kall + randn(N,N,Nreps) + i*randn(N,N,Nreps);  % Complex Noise, sigma=1

% -- 1. Convert to "reference" image and measure SNR
imall = 0*kall;
for n=1:Nreps
  imall(:,:,n) = (1/N)*ft(kall(:,:,n));	% FT all images
end;

im1 = abs(imall(:,:,1));
roisig = mean(im1(roi));	% Get mean signal in foreground
bgmean = mean(im1(bgroi));	% Get mean background
snrest = roisig / (bgmean/sqrt(pi/2));	% Estimate SNR
tt=sprintf('Estimated SNR is %g',snrest); disp(tt);

figure(1);
[mn,sd] = lec7kimhist(kall(:,:,1),'Reference',abs(im1(roi)));
refsnr=mn/sd;

% -- 2. Repeat measurement with multiple reps

pixsigs = imall(N/2,N/2,:);	% Get one pixel, all Nreps signals

figure(2); 
[mn,sd] = lec7kimhist(mean(kall(:,:,:),3),'Replicas',abs(pixsigs),refsnr);


% 3. Reduce FOV and measure Image SNR
% NOte that we will reconsruct the the original FOV here, could crop in y.

krFOV = kall(:,:,1);		% One k-space here.
krFOV(2:2:end,:)=0;		% Discard even lines (set to 0)
imrFOV = (2/N)*ft(krFOV);    % Normalize for N x (N/2)	 

figure(3); 
[mn,sd] = lec7kimhist(krFOV,'reduced FOV',abs(imrFOV(roi)),refsnr);


% 4. Reduce resolution and measure Image SNR
krres = kall(:,:,1);		% One k-space here.
krres = zpadcrop(krres,[N/2,N/2]);% Crop k-space
kspr = zpadcrop(ksp,[N/2,N/2]);   % Crop denoised k-space for ROI

imrres = (1/N)*ft(krres);    % Same normalizing (same k-space)
imrroi=0*imrres;
imrroi(N/4-N/16+3:N/4+N/16-2, N/4-N/16+3:N/4+N/16-2) = ones(N/8-4,N/8-4)*siglevel;
rroi = find(abs(imrroi(:))>siglevel/2);

figure(4); 
[mn,sd] = lec7kimhist(krres,'reduced res',abs(imrres(rroi)),refsnr);
%subplot(2,2,4);
%imroicheck = 0*imrroi; imroicheck(rroi)=1; 
%dispim(imroicheck);

% 5. Do averaging and measure Image SNR
kavg = mean(kall(:,:,1:2),3);	% Average 2 k-spaces
imavg= (1/N)*ft(kavg);		% Normalize NxN

figure(5); 
[mn,sd] = lec7kimhist(kavg,'Averaged(2x)',abs(imavg(roi)),refsnr);



