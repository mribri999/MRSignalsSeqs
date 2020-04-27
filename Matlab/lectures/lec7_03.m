
% Lecture 7, Example 03
% 
% Sample central half of k-space twice, outer half once.
%

N = 256;
ksp = randn(1.5*N,N)+i*randn(1.5*N,N);	% complex, sigma=1
ksp(N/4+1:3*N/4,:) = 0.5* ( ksp(N/4+1:3*N/4,:) + ksp(N+1:end,:)); % Average
ksp = ksp(1:N,:);	% Crop

figure(1); 
set(0,'defaultAxesFontSize',20);        % Default font sizes
dispim(ksp);
title('k-space');

im = 1/N*ft(ksp);
figure(2); 
set(0,'defaultAxesFontSize',20);        % Default font sizes
dispim(im);
title('Image Noise');

figure(3);
mn = mean(real(im(:)));
sd = std(real(im(:)));
ghist(real(im(:)),0,1/sqrt(1.5),[],'Histogram of Real(noise)','Data','\sigma = 1/sqrt{1.5}');
tt=sprintf('Mean and Std.Dev are %g and %g', mn, sd); disp(tt);
tt=sprintf('Std.Dev * sqrt(1.5) is %g',sd*sqrt(1.5)); disp(tt);
tt=sprintf('sqrt(9/8) = %g',sqrt(9/8)); disp(tt);



