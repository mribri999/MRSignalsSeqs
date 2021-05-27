%% This script demonstrates a random walk in one and three dimensions.
%
% Rad229_Lecture15_RandomWalk.m
%
% DBE@STANFORD.EDU (May 2020) for Rad229
%
% Example demos:
% 1) Examine the random walk stepping code.
%
% 2) Change the temperature.
%
% 3) Introduce anisotropy.
%
% To Do:
%  Define D and pick dx etc., or pick dx and compute D


%% This function demonstrates a 1D random walk
N_step = 100;               % Number of random walk steps
x = zeros( 1 , N_step);     % Starting x-position
dx = 55.5e-6;               % Spatial step size per dt [mm]
dt = 1e-6;                  % Temporal step size [s]

rng('default');             % Set the random-number-generator seed
p_xp = 0.5;                 % Probability (P) of a +x-step
% p_xm = 1 - p_xp;          % Probability (P) of a -x-step

for n = 2:N_step            % Loop over number of steps
  step = 2*(rand>p_xp)-1;   % +/- one step
  x(n) = x(n-1) + dx*step;  % Add +/- displacement step
end
t=0:dt:dt*(numel(x)-1);

figure; plot(t,x); title('1D Random Walk');
  xlabel('Time [s]'); ylabel('X-position [mm]');
  set(get(gca,'XAxis'),'Exponent',-6);
  set(get(gca,'YAxis'),'Exponent',-6);
  
%% This function demonstrates a 3D random walk
N_spin = 49;          % Number of spins to simulate
temp = 37;            % Sample temperature [C]

% Time-scale parameters
dt = 25e-6;                % Time step [s]
T = 10e-3;                 % Simulation duration [s]
N_step = round(T/dt);      % Number of time steps

% Length scale parameters
D_H20 = Diffusion_Coeff_H20; % [mm^2/s]
dx = sqrt(2*D_H20(temp)*dt); % Step size [mm] depends on diffusion coeff. and temperature
dy = 2*dx;                     % Isotropic when dy = dx
dz = 4*dx;                     % Isotropic when dz = dx
dr = cat(3,dx.*ones(N_spin,N_step),dy.*ones(N_spin,N_step),dz.*ones(N_spin,N_step));

% Define a uniform random step direction (sign)
tmp=rand(N_spin,N_step,3);   % Uniform random array (space x time x 3D)
SIGN=2*(tmp>0.5)-1;          % Unit (+/-) steps

% Define a uniform random step in x, y, OR z
step=rand(N_spin,N_step,3);            % Uniform random array of steps(space x time x 3D)
dR=SIGN.*dr.*(step==max(step,[],3));   % Step *only* in the direction of the largest random number

% Accumulate the individual steps into a path
path=cumsum(dR,2);

% Make a figure -OR- a movie
set(groot,'defaultLineLineWidth',2,'defaultLineMarkerSize',25,'defaultAxesFontSize',18);
f=figure; hold on; axis equal square; rot3d;
  set(f,'Units','Centimeters','Position',[10 10 14 12]);
title(['Diffusion Random Walk @',num2str(temp,'%02d'),'C']);
xlabel('X-position [mm]'); lab(1)=ylabel('Y-position [mm]');
axis([-0.01 0.01 -0.01 0.01 -0.01 0.01]);
path_color = hsv(N_spin);

% % Make a FIGURE
for k=1:N_spin
  p(k)=plot3(path(k,:,1),path(k,:,2),path(k,:,3));
  set(p(k),'color',path_color(k,:));
end
% saveas(f,['./Figures/RandomWalk_',num2str(temp,'%02d'),'_Final'],'png');

% % Make a MOVIE
% for t=1:N_step
%   for k=1:N_spin
%     p(k)=plot3(path(k,1:t,1),path(k,1:t,2),path(k,1:t,3));
%     set(p(k),'color',path_color(k,:));
%   end
%   saveas(f,['./Figures/RandomWalk_',num2str(temp,'%02d'),'_',num2str(t,'%03d')],'png');
%   delete(p);
% end

% % Plot a histograms...
% figure; hist(path(:,end,1),[min(path(:,end,1)):5:max(path(:,end,1))]);
% figure; hist(path(:,end,2),[min(path(:,end,2)):5:max(path(:,end,2))]);
% figure; hist(path(:,end,3),[min(path(:,end,3)):5:max(path(:,end,3))]);