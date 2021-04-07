% Script to generate some sample saturation curves for MT

clc;
clear;
close all;

TR = 40;
alpha = 15;
phi0 = 117*pi/180;
T1_MT = [1200 1200];
f_MT = 0.35;
k_MT = 4.3e-3;
T2_MT = 350;
T2b=12e-6;
[ff,G] = SuperLorentzian(T2b);
f_offsets = [0.5 [1:1:25]]*1E3;
b1s = [15:5:35];

for zz=1:length(b1s)
b1 = b1s(zz);
gam = 267.5221 *1e-3; % rad /ms /uT
trf = (pi*alpha/180)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;

npulse=300;
phi = RF_phase_cycle(npulse,phi0);

[s0,P] = epg_X_rfspoil(alpha*pi/180,phi,T1_MT,T2_MT,TR,f_MT,k_MT,pi*gam^2*b1sqrdtau*15.1*1E-3);



for ii=1:length(f_offsets)
    f_off = f_offsets(ii);
    GG = interp1(ff,G,f_off);
    WT = pi*gam^2*b1sqrdtau*(15.1+GG)*1E-3; 
    [s,P] = epg_X_rfspoil(alpha*pi/180,phi,T1_MT,T2_MT,TR,f_MT,k_MT,WT);
    s_all(:,ii) = abs(s);
    
end
if zz==1
figure
subplot(1,2,1)
selected=1:5:25;
plot(abs(s0))
xlabel('TR number')
hold on
plot(s_all(:,selected),'linewidth',1.5);
legends={};
for ii=1:length(selected)
    legends{ii+1}=['offset frequency=' num2str(f_offsets(ii))];
end
legends{1} = 'S0';
legend(legends)
end

subplot(1,2,2)
amp_s0 = abs(s0(end));
amp_sall = s_all(end,:);
plot(f_offsets/1000,amp_sall/amp_s0,'linewidth',1.5)
xlabel('Frequency offset[KHz]')
ylabel('S/S_0')
hold on
end
legend('15uT','20uT','25uT','30uT','35uT')