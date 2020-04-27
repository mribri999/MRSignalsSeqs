
% Lecture 7, Example 01
% 
% Calculate noise in k-space, image space, then take magnitudes
%
nsig = 1;	% Noise sigma parameter (real and imaginary)
N = 256;	% Image/k-space size.

kr = randn([N,N])*nsig;		% Generate gaussian noise (real part)
ki = randn([N,N])*nsig;		% Generate gaussian noise (imag part, same)
k = kr+i*ki;			% Combine

tt=sprintf('Std.Dev of noise in k-space is %g (real) and %g (imag)', ...
	std(kr(:)),std(ki(:)));
disp(tt);		% Display calculated stats

im = N * ift(k);	% iFFT with scaling of sqrt(N*N)=N

tt=sprintf('Std.Dev of noise in k-space is %g (real) and %g (imag)', ...
	std(real(im(:))),std(imag(im(:))));
disp(tt);

% -- Calculate noise as a function of SNR, with magnitude images

s=[0:0.1:10];		% s = Signal, so same as SNR if nsig==1
for p=1:length(s)
  mn(p) = mean(abs(im(:)+s(p)));	% Mean of magnitude signal
  sd(p) = std(abs(im(:)+s(p)));		% Std.deviation of magnitude signal
end;

set(0,'defaultAxesFontSize',20);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width

figure(1);			
dispim(cat(2,k,im)); axis equal;	% Display
title('k-space and image noise');	% 
axis off;

figure(2);
plot(s,mn-s,'b--',s,sd,'r-'); 		% Plot mean ans std.deviation
hold on;
plot(0,sqrt(pi/2)*nsig,'kx');		% Plot expectd mean
plot(0,sqrt(2-pi/2)*nsig,'ko');		% Plot expected std.dev
hold off;
leg1=legend('Mean Bias','Std.Deviation','sqrt(pi/2)','sqrt(2-pi/2)');
lplot('Signal','Statistics','Magnitude Noise Statistics vs Signal');
setprops;



