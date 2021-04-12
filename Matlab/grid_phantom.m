% This function creates a simple grid-pattern phantom.
%
% SYNTAX - grid_im = grid_phantom(imsize)
%
% INPUT - imsize - integer that defines square matrix size
%
% OUTPUT - grid_im - Image of a circular phantom with a grid pattern
% 
% Original code creation by Michael Loecher. DBE@STANFORD.EDU (April 2020 for Rad229)

function grid_im = grid_phantom(imsize)

if nargin == 0, imsize = 256; end

phantom_rad = 0.46;
phantom_mag = 0.9;
comb_os = 20;
comb_cutoff = -.8; % Makes initial comb function, between -1 and 1, smaller is thinner lines
blur_mod = 5; 
grid_freq = 20; % How many lines in the total image (this will get cropped)
grid_crop = 0.5; % Crop radius of grid (>0.5 = no crop)

im_background = zeros(imsize);
[XX, YY] = meshgrid( 1 : imsize , 1 : imsize);
XX = (XX - imsize / 2.0) / imsize;
YY = (YY - imsize / 2.0) / imsize;
RR = sqrt(XX.^2 + YY.^2);

im_background(RR < phantom_rad) = phantom_mag;

Ncomb = imsize*comb_os;
cos_comb = cos(2.0*pi.*grid_freq.*linspace(-0.5, 0.5, Ncomb));
comb = ones(1, Ncomb);
comb(cos_comb < comb_cutoff) = 0.0;

comb = smoothdata(comb,'gaussian',blur_mod*comb_os);
comb = comb(1:comb_os:end);
comb = comb'*comb;

grid_mod = comb;
grid_mod(abs(XX) > grid_crop) = 1.0;
grid_mod(abs(YY) > grid_crop) = 1.0;

grid_im = im_background.*grid_mod;

end

