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
ksamp = ones(256,256);  % FULLY sampled

% == Parallel Imaging with Calibration
if (0)
    ksamp(1:2:end,:) = 0;   % <<=== MODIFY TRAJECTORY HERE!!
    kcalhw = 0;                        % Calibration half-width
    ksamp(129-kcalhw:128+kcalhw,:)= 1;    % Add in calibration region
end

% == Random Sampling
if (1) 
    krand = rand(256,256);  % Random numbers 0 to 1
    R=4;
    ksamp = 0*ksamp;        % Default is not sampled.
    ksamp(find(krand(:)<(1/R)))=1;   % Choose sample locations
    kcalhw = 4;                        % Calibration half-width
    ksamp(129-kcalhw:128+kcalhw,129-kcalhw:128+kcalhw)= 1;    % Add in calibration region
end
kactual = kdata.*ksamp; % 

im = psf2d(ksamp);    % Default, 256x256 k-space.

figure;
dispim(ft(kactual));  
title('Image');






