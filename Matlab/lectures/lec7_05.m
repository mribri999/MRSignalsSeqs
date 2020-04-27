
% Lecture 7, Example 05
% 
% A simple 2-coil example of optimal weights.

S = [3;1];	% Sensitivities of 2 coils (at a pixel)

% -- Test all relative weights of these coils

w1 = [0:0.01:0.95];
w2 = 1-w1;

for k = 1:length(w1)
  W = [w1(k) 1-w1(k)];
  sig = W*S;		% Signal
  nvar = sum(W.^2);	% Variance
  snr(k) = sig/sqrt(nvar);	%
end;


[y,yi] = max(snr);
plot(w1,snr,'r-',w1,w1./w2,'b--',w1(yi)*[1 1],[0 4],'k-');
lplot('Weight on Coil 1','SNR','SNR vs weights');
legend('SNR','W1/W2','optimal');
axis([0 1 0 1.25*y]);

[W,g,imn] = senseweights(S);
disp(['SENSE Weights']); disp(W);
SENSEsnr = (W*S)/sqrt(imn);
disp('SENSE SNR'); disp(SENSEsnr);





