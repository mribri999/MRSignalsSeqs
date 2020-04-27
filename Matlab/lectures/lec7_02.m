
% Lecture 7, Example 02
% 
% Plot noise distributions.
%

im = randn(256,256)+i*randn(256,256);	% complex, sigma=1

figure(1);
ghist(real(im(:)),0,1,[],'Histogram of Real(noise)')

figure(2);
ghist(abs(im(:)),[],[],[],'Histogram of abs(noise)')

figure(3);
ghist(abs(im(:)+4),[],[],[],'Histogram of abs(4+noise)')

