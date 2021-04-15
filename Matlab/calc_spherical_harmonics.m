% This function computes the the spherical harmonic coefficients that are
% specific to a scanner and used to compute the non-linear gradient induced
% displacement distortion.
%
% SYNTAX - B = calc_spherical_harmonics(Alpha, Beta, X, Y, Z, R0)
%
% INPUTS - Alpha - Coefficient for the cosine terms
%          Beta - Coefficient for the sine terms
%          X - X-position at which to compute the harmonic [1 x m]
%          Y - Y-position at which to compute the harmonic [1 x m]
%          Z - Z-position at which to compute the harmonic [1 x m]
%          R0 - Scaling factor
%
% OUTPUTS - B - Spherical harmonic values at (X,Y,Z) [1 x m]
%
% Original code creation by Michael Loecher. DBE@STANFORD.EDU (April 2020 for Rad229)

function B = calc_spherical_harmonics(Alpha, Beta, X, Y, Z, R0)
% Calculate displacement field from spherical harmonic coefficients

nmax = size( Alpha , 1 ) - 1;

% Convert to spherical coordinates
R = sqrt( X.^2 + Y.^2 + Z.^2);
R_eps = R + .0001;
Theta = acos( Z ./ R_eps );
Phi = atan2( Y ./ R_eps , X ./ R_eps );

% Evaluate the Legendre polynomial with normalization
B = 0;
for n = 0 : nmax
  P = gradient_nonlin_legendre(n,cos(Theta));
  F = ( R / R0 ) .^ n;
  for m = 0 : n
    F2 = Alpha( n + 1 , m + 1 ) * cos( m * Phi ) + Beta( n + 1 , m + 1) * sin( m * Phi );
    B = B + F .* P( m + 1 , : ) .* F2;
  end
end

return

function P = gradient_nonlin_legendre( n , X )
% Calculate the Legendre polynomial with normalization

P = legendre(n,X);

for m = 1 : n
  normfact = (-1) ^ m * sqrt( (2 * n + 1) * factorial( n - m ) / ( 2 * factorial(n + m) ) );
  P( m + 1 , : ) = normfact * P( m + 1 , : );
end

return