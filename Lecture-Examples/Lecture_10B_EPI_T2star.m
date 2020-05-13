%% This script demonstrates T2* blurring echo planar imaging (EPI).
%
% DBE@STANFORD.EDU (May 2020) for Rad229

%% Default parameters
N=128;                    % Matrix is NxM (then padded later to accomodate motion)
T2star = [10e-3 25e-3 50e-3 100e-3 1000e-3]; % Range of T2-star values
ESP = 1e-3; % EPI echo spacing

%% Create a resolution phantom with some noise
Obj0 = zeros(N);
C0 = 16:2:N;  C1 = 16:2:N; R0 = 9:24:N;  R1 = 25:24:N;
for c=1:5, for r=1:5
    Obj0(R0(r):R1(r),C0(c):C1(c))=T2star(r); 
end; end

C0 = 48:4:N;  C1 = 49:4:N; R0 = 9:24:N;  R1 = 25:24:N;
for c=1:5, for r=1:5
    Obj0(R0(r):R1(r),C0(c):C1(c))=T2star(r); 
end; end

C0 = 88:8:N;  C1 = 91:8:N; R0 = 9:24:N;  R1 = 25:24:N;
for c=1:5, for r=1:5
  Obj0(R0(r):R1(r),C0(c):C1(c))=T2star(r); 
end; end

Obj0=Obj0'; % Want the resolution elements parallel to phase encode direction

%% Create a time series of phantoms with T2-star decay and noise
Obj = repmat(Obj0,[1 1 N]);
for n=1:N
  Obj(:,:,n)=exp(-(n*1e-3)./Obj0);
end

%% Create a noisy object and k-space
Obj = Obj + (randn(size(Obj))+1i*randn(size(Obj)))/50; % Add some complex noise
kspc = zeros(size(Obj));
for n=1:N
  kspc(:,:,n)=fftshift(fft2(Obj(:,:,n))); % Define the k-space
end

%% Simulate the TOP-DOWN EPI trajectory...grab one k-line from each time, sign reverse
k_epi = zeros(N);
for n=1:N
  k_epi(n,:) = kspc(n,:,n); % Grab one echo per ky-line during T2-star decay
end

%% Simulate a CENTER-OUT EPI trajectory
k_cntout = zeros(N);
n_cntout = lin([(N/2+1):1:N; N/2:-1:1]);
ind=1;
for n=1:N
  k_cntout(n_cntout(n),:) = kspc(n_cntout(n),:,ind); % Grab one echo per ky-line during T2-star decay
  ind=ind+1;
end

%% Sign-reversed for the odd echoes
% for k=1:2:N, k_epi(k,:)=fliplr(k_epi(k,:)); end

IM_epi = ifft2(fftshift(k_epi)); % Image of object
IM_cntout = ifft2(fftshift(k_cntout)); % Image of object

%% Create an EPI trajectory figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0)); axis equal tight;
 title('Object');        caxis([0 0.075]);
subplot(2,2,2); imagesc(abs(IM_epi)); axis equal tight;
  title('Image');         caxis([0 0.5]);
subplot(2,2,3); imagesc(abs(kspc(:,:,1))); axis equal tight;
  title('k-space (Obj)'); caxis([0 150]);
subplot(2,2,4); imagesc(abs(k_epi)); axis equal tight;
  title('k-space (EPI)'); caxis([0 300]);
  
%% Create an Center-Out EPI trajectory figure...
figure; hold on; colormap(hot);
subplot(2,2,1); imagesc(abs(Obj0)); axis equal tight;
 title('Object');        caxis([0 0.075]);
subplot(2,2,2); imagesc(abs(IM_cntout)); axis equal tight;
  title('Image');         caxis([0 0.5]);
subplot(2,2,3); imagesc(abs(kspc(:,:,1))); axis equal tight;
  title('k-space (Obj)'); caxis([0 150]);
subplot(2,2,4); imagesc(abs(k_cntout)); axis equal tight;
  title('k-space (EPI)'); caxis([0 300]);
