% Lecture 13, Example 03
%
% Homodyne reconstruction

set(0,'defaultAxesFontSize',14);        % Default font sizes
figure(1);

% -- Generate Sample Data (Array of squares)
[kx,ky] = meshgrid([-128:127]/256,[-128:127]/256);
[cx,cy] = meshgrid(2*[-2:2],2*[-2:2]);
centers = [cx(:)+i*cy(:); [-8; -6; 6; 8]; i*[ -8; -6; 6; 8]]*10;
kdata = ksquare(centers,19,kx+i*ky,0,0); % k-space data 

im = ft(kdata);	% -- Reconstruct the image

% -- Run the Homodyne example code
imhd = homodyneshow(im,16);

