% -- Show distributions Gaussian, Rician (& Rayleigh)
%

N=100;					% # histogram bins
Ns=500;					% # samples
if (~exist('sigmean'))
  sigmean = 0;				% Mean signal
end;

dat = randn(Ns,Ns)+i*randn(Ns,Ns)+sigmean;	
dat = dat(:);
[yreal,x] = hist(real(dat),N);
[ymag]= hist(abs(dat),x);
mn = mean(abs(dat));
sg = std(abs(dat));
hold off;
bar(x,[ymag' yreal'],1);
hold on;
tt = sprintf('%g*normpdf(x,%g,%g)',sqrt(2*pi*sg*sg)*max(ymag),mn,sg);
fplot(tt,[min(x),max(x)],':');

legend('Rician/Rayleigh (|sig|)','Gaussian (Real{sig})','Gaussian w/ Rician stats');
grid on;
title('Signal Distributions, mean=0');
xlabel('Signal');
ylabel('Frequency');

setprops;

