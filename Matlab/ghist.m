%
%	Function ghist(data,gmean,gsig,nbins,gtitle,gleg1,gleg2)
%
%	Function makes is histogram of a signal, and superimposes
%	a gaussian fit.
%
%	INPUT:
%		data = signal to make histogram for.
%		gmean = estimate for gaussian fit, mean [default=mean(data)]
%		gsig = estimate for gaussian fit, sigma [default=std(data)]
%		nbins = number of bins
%		gtitle = plot title
%		gleg1 = legend top line (describe data)
%		gleg2 = legend second line (describe gaussian)
%
function ghist(data,gmean,gsig,bins,gtitle,gleg1,gleg2)

N = length(data(:));

if nargin < 2 || length(gmean)==0 gmean = mean(data); end;
if nargin < 3 || length(gsig)==0 gsig = std(data); end;
if nargin < 4 || length(bins)==0 
	nbins=100; 
	bins = 4*gsig*[-1:2/(nbins-1):1];	% Histogram bins
end;
if nargin < 5 || length(gtitle)==0 gtitle='Data and Gaussian'; end;
if nargin < 6 || length(gleg1)==0 gleg1 = 'Data'; end;
if nargin < 7 || length(gleg2)==0 gleg2 = 'Gaussian Fit'; end;



hold off;
histogram(data(:),bins); 	          	% Plot Histogram lplot('Signal','Likelihood',gtitle,[min(bins),max(bins),0,10*sqrt(N)]);
hold on;
gscale = N/length(bins)/1.5;
%gscale = N/length(bins)/sqrt(sqrt(bsigma));             % Scaling
plot(bins,gscale*gaussian(bins,gmean,gsig),'r-');  % Plot normal distribution.
legend(gleg1,gleg2);
hold off;
grid on;




