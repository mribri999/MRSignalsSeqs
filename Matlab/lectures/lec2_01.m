
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
dispim(cat(2,im1,im2,im3)); 
title('Images for a, b and c');
axis equal;
axis off;

