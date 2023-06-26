% Script to compare to original implementation from
% https://github.com/mriphysics/EPG-X. Need to add such library to the
% Matlab path.

clc;
clear;
close all;

TR = 5;
alpha = 10;
phi0 = 117*pi/180;
T1_MT = [779 779];
f_MT = 0.117;
k_MT = 4.3e-3;
T2_MT = 45;

%%% Relaxation parameters: exchange
T1x = [1000 500];
T2x = [100 20];    
kx = 2e-3;
fx = 0.2;


G = 15.1;% us 
b1 = 13; % uT
gam = 267.5221 *1e-3; %< rad /ms /uT
trf = (pi*alpha/180)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;

npulse=250;
phi = RF_phase_cycle(npulse,phi0); %function taken from Malik et al.
WT = pi*gam^2*b1sqrdtau*G*1e-3; 

[s,P] = epg_X_rfspoil(alpha*pi/180,phi,T1_MT,T2_MT,TR,f_MT,k_MT,WT);
[s2,P] = epg_X_rfspoil(alpha*pi/180,phi,T1_MT,T2_MT,TR,f_MT*2,k_MT,WT);
[s3,P] = epg_X_rfspoil(alpha*pi/180,phi,T1_MT,T2_MT,TR,f_MT*3,k_MT,WT);


figure
subplot(2,1,1)
plot(abs(s(:,:).'),'linewidth',2)
hold on
plot(abs(s2(:,:).'),'linewidth',2)
plot(abs(s3(:,:).'),'linewidth',2)
legend('f=0.12','f=0.23','f=0.35')
xlabel('TR number')

title('Own implementation')
ylim([0 0.17])

%% original
[smt,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),TR,T1_MT,T2_MT,f_MT,k_MT,G); %function taken from Malik et al.
[smt2,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),TR,T1_MT,T2_MT,f_MT*2,k_MT,G);
[smt3,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),TR,T1_MT,T2_MT,f_MT*3,k_MT,G);


subplot(2,1,2)
plot(abs(smt),'linewidth',2)
hold on
plot(abs(smt2),'linewidth',2)
plot(abs(smt3),'linewidth',2)
title('Original')
ylim([0 0.17])
legend('f=0.12','f=0.23','f=0.35')
xlabel('TR number')