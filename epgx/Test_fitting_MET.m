%% File to simulate T2 fitting with different exchange rates

clc;
clear;
close all;


run('/home/gustavo/Downloads/qMRLab-master/startup.m') %download from https://qmrlab.org/
warning ('off','all');

T1 = [1000 500];
T2 = [120 20];    
list_kxs = [0:1:20]*1e-3; % list of exchange rates to test;
list_f = [1:40]*1e-2; % list of fractions to test;
 
for jj=1:length(list_kxs)
    ka = list_kxs(jj);
    for ii=1:length(list_f)
        f = list_f(ii);
        esp=5;
        angles = pi*ones(1,50);
        Npulse = length(angles);
        [s,P] = epg_X_CMPG(angles,f,T1,T2,esp,ka,0);
        all_s(ii,jj,:) = abs(s(3,:));
    end
end

% load presaved model to use for fitting
Model = qMRloadObj('mwf.qmrlab.mat');
Model.Prot.MET2data.Mat = esp*(1:Npulse)';
data = struct();
data.MET2data= reshape(all_s,[length(list_f) length(list_kxs) 1 Npulse]);
data.Mask= ones([length(list_f) length(list_kxs) 1]);
FitResults = FitData(data,Model,0);
save('FitResultsMET')

figure
subplot(1,3,1)
selected = 1:2:length(list_kxs);
plot(list_f,list_f);
hold on
plot(list_f,FitResults.MWF(:,selected)/100,'*','linewidth',1)
xlabel('real fraction')
ylabel('estimated fraction')
legends={};
for ii=1:length(selected)
    legends{ii+1}=['ka=' num2str(list_kxs(selected(ii)))];
    
end
legends{1} = 'Ground truth';
legend(legends);
xlim([0 0.4])
ylim([0 0.4])
title('Fraction estimation')

subplot(1,3,3)
selected = 1:2:length(list_kxs);
plot(list_f,T2(2)*ones(size(list_f)));
hold on
plot(list_f,FitResults.T2MW(:,selected),'*','linewidth',1)
xlabel('fraction')
ylabel('Estimated T2_b (ms)')
legends={};
for ii=1:length(selected)
    legends{ii+1}=['ka=' num2str(list_kxs(selected(ii)))];
    
end
legends{1} = 'Ground truth';
title('T2_b estimation')

subplot(1,3,2)
selected = 1:2:length(list_kxs);
plot(list_f,T2(1)*ones(size(list_f)));
hold on
plot(list_f,FitResults.T2IEW(:,selected),'*','linewidth',1)
xlabel('fraction')
ylabel('Estimated T2_a(ms)')
legends={};
for ii=1:length(selected)
    legends{ii+1}=['ka=' num2str(list_kxs(selected(ii)))];
    
end
legends{1} = 'Ground truth';
% legend(legends);
title('T2_a estimation')
