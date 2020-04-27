
% Lecture 7, Example 07
% 
% Estimate noise covariance matrix.
%
psi = [1 .3 .1; .3 2 .2; .1 .2 1.5];
N = 5000;

n = mvnrnd(zeros(1,3),psi,N)+i*mvnrnd(zeros(1,3),psi,N);

psi_est = 0;
for k=1:N
  psi_est = psi_est + n(k,:)'*(n(k,:));
end;
 
psi_est = 1/2/N * psi_est

[psi psi_est]


