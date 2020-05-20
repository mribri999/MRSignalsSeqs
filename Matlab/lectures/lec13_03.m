% Lecture 13, Example 03
%
% Plot point-spread functions using psf2d.m

set(0,'defaultAxesFontSize',14);        % Default font sizes
figure(1);


% Define Some Goals here...!
disp('TRICKS - annular sub-sampling circ.m');
disp('Random sampling');
disp('Parallel Imaging with a calibration region');
disp('TSENSE Sampling - odd ky then even, look at psd2d time option');


% -- Generate Sample Data (Array of squares)
[kx,ky] = meshgrid([-128:127]/256,[-128:127]/256);
[cx,cy] = meshgrid(2*[-2:2],2*[-2:2]);
centers = [cx(:)+i*cy(:); [-8; -6; 6; 8]; i*[ -8; -6; 6; 8]]*10;
kdata = ksquare(centers,19,kx+i*ky,0,0); % k-space data 
% -- Could use gridding to add to a grid, but better way is below...
%kmatrix = gridmat(kx+i*ky,kdata,ones(size(kdata)),256);

% ==== DEFINE SAMPLING TRAJECTORY HERE ====
ksamp = ones(256,256);	% FULLY sampled

ksamp(1:2:end,:) = 0;	% <<=== MODIFY TRAJECTORY HERE!!

kactual = kdata.*ksamp;	% 

im = psf2d(kactual);	% Default, 256x256 k-space.

