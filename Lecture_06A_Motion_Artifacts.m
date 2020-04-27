%% This script demonstrates the impact of steady and pulsatile motion on
%  both images and k-space.
%
% DBE@STANFORD.EDU (April 2020) for Rad229

%% Set-up some default parameters
N=128;               % Matrix is NxN (then padded later to accomodate motion)
dt=10e-3;            % Time step
N_step=16;           % Number of time steps
dur=dt*(N_step-1);   % Duration of motion
t=0:dt:dur;          % Time vector during motion

%% Pulsatile "vessel" signal frequency and amplitude
% Small vessel; non-pulsatile
f_p=0;     % Pulse frequency
A_p=0;     % Pulse amplitude
pix_p=1;   % Pixel "radius"

% Intermediate vessel; pulsatile
% f_p=4;   % Pulse frequency
% A_p=10;  % Pulse amplitude
% pix_p=3; % Pixel "radius"

%% Define the motion field
dy_max=round(0);     % No Y-displacement
% dy_max=round(6);     % Small y-displacement
% dy_max=round(48);    % Large y-displacement
yt=dy_max.*([0.00 0.00 0.00 0.00 0.33 0.33 0.33 0.33 0.67 0.67 0.67 0.67 1.00 1.00 1.00 1.00]);

% No X-displacement
dx_max=round(0);     % Max x-displacement
xt=dx_max.*([0.00 0.00 0.00 0.00 0.33 0.33 0.33 0.33 0.67 0.67 0.67 0.67 1.00 1.00 1.00 1.00]);

%% Generate the phantom and apply the motion
Obj = phantom('Modified Shepp-Logan',N);
P = padarray(zeros(N,N,N_step),[abs(dy_max) abs(dx_max)],0,'both');  % Pad for motion

for n=1:length(t)
  p_ind=(N/2-pix_p):(N/2+pix_p);
  pulse_amp(n)=A_p*(1+sin(2*pi*f_p*(n-1)/N_step));
  Obj(p_ind,p_ind)=pulse_amp(n);

  P(1+yt(n):N+yt(n),1+xt(n):N+xt(n),n)=Obj;
  P(:,:,n)=P(:,:,n)+(randn(size(P(:,:,n)))+i*randn(size(P(:,:,n))))/25; % Add some complex noise
end

%% Define k-space for each phase
kspc = complex(zeros(size(P)));
for n=1:length(t)
  kspc(:,:,n)=fftshift(fft2(P(:,:,n)));
end

%% Extract a few ky-lines from each time step during motion
VPS=size(P,1)/N_step;  % Views-per-segment (lines of k-space per time step)
for m=1:N_step
  ind0=1+(m-1)*VPS;  % Index to first ky-line
  ind1=m*VPS;        % Index to last ky-line
  kspc_motion(ind0:ind1,:)=kspc(ind0:ind1,:,m); % Copy lines from the mth phase
  size(kspc_motion)
end

Obj_blur=ifft2(fftshift(kspc_motion)); % Image of moving object

%% Generate figures...
set(groot,'defaultLineLineWidth',3,'defaultLineMarkerSize',25,'defaultAxesFontSize',15);
figure; hold on; colormap(hot);
subplot(2,2,1); hold on;
  plot(t,xt); plot(t,yt,'--'); plot(t,pulse_amp);
  title('Motion Path'); legend('X-path','Y-Path','Pulse');

subplot(2,2,2); hold on;
  imagesc(abs(kspc_motion)); caxis([0 max(abs(kspc_motion(:)))/25]);
  title('Motion k-space MAG'); axis equal tight;

subplot(2,2,3); hold on;
  imagesc(abs(P(:,:,1))); caxis([0 1]);
  title('Object Mag.'); axis equal tight;
  
subplot(2,2,4); hold on;
  imagesc(abs(Obj_blur)); caxis([0 1]);
  title('Motion Object'); axis equal tight;