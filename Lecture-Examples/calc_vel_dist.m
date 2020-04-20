%%
clearvars
clc
close all

%%
load('pcmri_data0.mat')
load('synth_coeff.mat');

im0 = mag_d;

[Vx, Vy, Vz] = ndgrid(1:size(im0,1), 1:size(im0,2), 1:size(im0,3));
[Sx, Sy, Sz] = transform_coords(Vx, Vy, Vz, coords);

[Dx, Dy, Dz] = get_distortion(Sx, Sy, Sz, Alpha_x, Alpha_y, Alpha_z, ...
                                         Beta_x, Beta_y, Beta_z, R0);

Sx_d = Sx + Dx;
Sy_d = Sy + Dy;
Sz_d = Sz + Dz;

[Vx_d, Vy_d, Vz_d] = transform_coords(Sx_d, Sy_d, Sz_d, inv(coords));
