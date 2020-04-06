
% Lecture 1, Example 01
% 
% Simple image reconstruction in matlab
%
load rawkneedata
%% a) Basic reconstruction
im1=ft(dat);

%% b) Skip every 2nd line
dat2=dat; dat2(:,1:2:end)=0;  % Zero out every 2nd line 
im2=ft(dat2);

%% c) Zero-out outer half of k-space
dat3=dat; dat3(:,1:64)=0; dat3(:,193:end)=0;   % Zero outer k-space
im3=ft(dat3);

% -- Display all images
figure(1);
dispim(im1); title('Basic Image'); axis off;

figure(2); 
dispim(im2);
title('Image with Every 2nd vertical (S/I) k-space Line'); axis off;

figure(3);
dispim(im3);
title('Image with Outer k-space set to 0'); axis off;


