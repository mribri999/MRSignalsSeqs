% Script to compare to original implementation from
% https://github.com/mriphysics/EPG-X. Need to add such library to the
% Matlab path.
clc;
clear;
close all;

f = 0.20; 
ka = 5e-2; 
T1 = [1000 500];
T2 = [100 20];
esp=5;
angles = pi*ones(1,50);
Npulse = length(angles);
[s,P] = epg_X_CMPG(angles,f,T1,T2,esp,ka,0);
figure
subplot(1,2,1)
plot(abs(s(:,:).'),'linewidth',2)
title('Own implementation')
legend('1st compartment','2nd compartment','voxel average','fontsize',13)

a0 = [90 180*ones(1,50)]*pi/180;

[s_gt,Fn_gt] = EPGX_TSE_BM(a0,esp,T1,T2,f,ka,'delta',0);  %function taken from Malik et al.
subplot(1,2,2)
plot(abs(s_gt))
plot(abs([squeeze(Fn_gt(52,:,:)) s_gt.']),'linewidth',2)
title('Original')
legend('1st compartment','2nd compartment','voxel average','fontsize',13)