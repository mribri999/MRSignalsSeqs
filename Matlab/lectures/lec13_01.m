% Lecture 13, Example 01
%
% Plot point-spread functions using psf2d.m

set(0,'defaultAxesFontSize',14);        % Default font sizes
figure(1);
im = psf2d;	% Default, 256x256 k-space.
subplot(2,2,1); title('k-space - Fully sampled');

figure(2);
k = ones(256,256);
k(1:3:end,:)=0;
k(2:3:end,:)=0;
im = psf2d(k); 
subplot(2,2,1); title('k-space - 3x undersampled ky');

figure(3);
ff = 1./(exp((abs([1:256]-128)-112)/16)+1);
ff = diag(ff); 
k = ones(256,256);
kf = ff*k*ff;	% Apply Fermi filter
im = psf2d(kf); 
subplot(2,2,1); title('k-space - Fermi-filtered');
