%**************************************************************************%
function B = calc_harmonics(Alpha, Beta, X, Y, Z, R0)
%MRIS_GRADIENT_NONLIN__SIEMENS_B  calculate displacement field from Siemens's coefficients

nmax = size(Alpha,1)-1;

% convert to spherical coordinates
R = sqrt(X.^2+Y.^2+Z.^2);
R_eps = R + .0001;
Theta = acos(Z./R_eps);
Phi = atan2(Y./R_eps,X./R_eps);


% evaluate the Legendre polynomial (using Siemens's normalization)
B = 0;
for n = 0:nmax
  P = mris_gradient_nonlin__siemens_legendre(n,cos(Theta));
  F = (R/R0).^n;
  for m = 0:n
    F2 = Alpha(n+1,m+1)*cos(m*Phi)+Beta(n+1,m+1)*sin(m*Phi);
    B = B+F.*P(m+1,:).*F2;
  end
end

end


%**************************************************************************%
function P = mris_gradient_nonlin__siemens_legendre(n,X)
%MRIS_GRADIENT_NONLIN__SIEMENS_LEGENDRE  normalize Legendre polynomial with Siemens's convention

P = legendre(n,X);

for m=1:n
normfact = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)));
P(m+1,:) = normfact*P(m+1,:);
end

 end