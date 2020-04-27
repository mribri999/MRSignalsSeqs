%
%	Function plots k-space, image and histogram, and returns mean and sd.

function [mn,sd] = lec7kimhist(ksp,ftitle,histpts,refsnr)

[N,M] = size(ksp);
im = (1/N)*ft(ksp);

subplot(2,2,1); dispim(log(1+abs(ksp))); 
tt = sprintf('kspace - %s',ftitle);
title(tt); axis off;
subplot(2,2,2); dispim(im); 
tt = sprintf('image - %s',ftitle);
title(tt); axis off;
subplot(2,2,3); 
[mn,sd] = ghist(histpts);
if (nargin < 4) refsnr = mn/sd; end;

title(ftitle);
tt=sprintf('%s:  mean=%g, stdev=%g, ratio=%g, relative=%g',ftitle,mn,sd,mn/sd,mn/sd/refsnr); disp(tt);

