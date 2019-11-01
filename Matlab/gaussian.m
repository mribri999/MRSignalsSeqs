%	Function y = gaussian(x,mn,sig)
%
%	Function generates the gaussian distribution output y 
%	based on x and the mean (mn) and sigma (sig).  x can be a 
%	number or array.
%
function y = gaussian(x,mn,sig)

y = exp(-(x-mn).^2 / (2*sig^2)) / (sqrt(2*pi)*sig);


