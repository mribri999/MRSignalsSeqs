function [Dx, Dy, Dz] = get_distortion(Sx, Sy, Sz, Alpha_x, Alpha_y, Alpha_z, Beta_x, Beta_y, Beta_z, R0)

x = Sx(:).';
y = Sy(:).';
z = Sz(:).';

bx = calc_spherical_harmonics(Alpha_x, Beta_x, x, y, z, R0);
by = calc_spherical_harmonics(Alpha_y, Beta_y, x, y, z, R0);
bz = calc_spherical_harmonics(Alpha_z, Beta_z, x, y, z, R0);

Bx = reshape(bx, size(Sx));
By = reshape(by, size(Sy));
Bz = reshape(bz, size(Sz));

Dx = Bx * R0;
Dy = By * R0;
Dz = Bz * R0;

end

