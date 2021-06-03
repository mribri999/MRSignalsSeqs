%% This script simulates the effect of diffusion encoding gradients on
%  spins that undergo a random walk.
%
% Rad229_Lecture15_RandomWalk_Diffusion.m
%
% DBE@STANFORD.EDU (May 2020) for Rad229
% Code base developed by Kévin Moulin
%
% ExG_maxe demos:
% 1) Examine the phase history amongst the spins. 
%    [Change the diffusion coefficient widely.]
%
% 2) Examine the effect of diffusion anisotropy.
%
% 3) Examine the effect of gradient direction.
%
% 4) EXTRA - delta, Delta, b-value changes...
%
% To Do:
%   Refine code overall
%   Check gradient waveform design
%   Compare phase variance to diffusion etc.
%   Add GRE ans SE options
%   Visualize the additional phase from diffusion (hidden by total phase)

%% Initialization
close all; clear all;
rng(0,'twister');                                                     % Set the random-number-generator seed
ANIM = 1;                                                           % Animation flag (slow...)

%% General parameters
sys = Rad229_MRI_sys_config;                                        % Define MRI system constants
World.DisplayPoint = 50;                                            % Number of simulated points to skip before display
World.Pixel = [0.5 0.5] .* 1e-3;                                      % Pixel dimension [m]

%% Diffusion parameters
diff.N_spins = 100;                                                 % Number of spins
diff.dim = 2;                                                       % Dimensions for spin diffusion [2D]
diff.D = [2.1 2.1].*1e-9;                                           % Coefficient of diffusion [(mm^2/s).*1e-6] = [m^2/s]
diff.dt = sys.dt;                                                   % Use the same dt as the global dt

%% Sequence parameters
seq.bvalue = 2000.* 1e6;                                            % b-value [(s/mm^2)*.1e6] = [s/m^2]
seq.delta = 12e-3;                                                  % Duration of one gradient lobe [s]
seq.Delta = 15e-3;                                                  % Duration between the gradient [s]
seq.direction = [1 1];                                              % Diffusion vector direction

Grad = generate_gradient_points_local( seq , World , sys );                   % List of gradient points with raster time corresponding to dt
Grad_dir = repmat( Grad , 1 , diff.dim ) .* seq.direction;                  % List of gradient points per physical axis

seq.duration = length(Grad);                                        % List of gradient points with raster time corresponding to dt

%% Init the Simulation
% Distribute N spins in the center part of the pixel for each dimension (For the 2D case, only X and Y)
spin = zeros( diff.N_spins , diff.dim , seq.duration );  % [N spins, N dimensions, N gradient raster points]
for cpt = 1 : diff.dim
  spin( : , cpt , 1 ) = 0.2 * randn(diff.N_spins,1) .* World.Pixel(cpt); % For X positions
end

% Initialize the MR signal per spin and per time point
MR_Phase_Molecule = zeros( diff.N_spins , seq.duration , 1);
MR_Signal_Pixel = complex( zeros( seq.duration , 1 ) , zeros( seq.duration , 1 ) );

%% Simulation
% Per time point we do the following:
% 1) Update the water molecule positions
% 2) Apply the gradients to the water molecules
% 3) Calculate the MR signal at the pixel scale
% 4) (Eventually) Update the plots

for cpt = 2 : 1 : seq.duration
  
  % 1) Move the spins according the diffusion parameters
  spin( : , : , cpt ) = spin( : , : , cpt-1 ) + diffusion_step_local( diff );
  
  % 2) Apply the gradients to the spin position and  accumulate the phase
  % The numerical integration is done with a cumulative sum (ALT: see trapezoid method)
  
  % 2a) Repeat the gradient direction vector per spin (convert the gradient to T/m)
  Grad_dir_tmp = repmat( Grad_dir( cpt , : ) , diff.N_spins , 1 ) * 1e-3;
  
  % 2b) Numerical integration of the spin positions x gradient, we sum over the dimensions (result in radian)
%   MR_Phase_Molecule(:,cpt)=MR_Phase_Molecule(:,cpt-1)+ nansum(World.Gamma *World.dt*Grad_dir_tmp.*squeeze(spin(:,:,cpt)),2);
  diff_phase = 2 * pi * sys.gamma_bar * sys.dt * Grad_dir_tmp;
  MR_Phase_Molecule( : , cpt ) = MR_Phase_Molecule( : , cpt-1 ) + sum( diff_phase .* squeeze( spin( : , : , cpt ) ) , 2 , 'OmitNaN' );
   
  % 3) Calculate the MR signal at the voxel scale
%   MR_Signal_Pixel=nanmean(cos(MR_Phase_Molecule(:,cpt))+1i*sin(MR_Phase_Molecule(:,cpt)));
  MR_Signal_Pixel = mean ( cos(MR_Phase_Molecule( : , cpt ) ) + 1i * sin( MR_Phase_Molecule( : , cpt ) ) , 'OmitNaN' );
  
end

%% 4) Display the simulation
figH = create_figure_local(MR_Phase_Molecule,spin,Grad_dir,seq.duration,World.Pixel,diff);
Grad_dir( end/2 : end , : ) = -Grad_dir( end/2 : end , : );

if ANIM
  %Uncomment to create animation/ comment to see only the final result
  for cpt = 2 : 1 : seq.duration
    if mod( cpt , World.DisplayPoint ) == 0
      update_figure_local( figH , MR_Phase_Molecule , spin , Grad_dir , cpt );
      pause(0.01);
    end
  end
end

update_figure_local(figH,MR_Phase_Molecule,spin,Grad_dir,cpt);

% Coefficient of diffusion measured

disp(['Diffusion attenuation = ' num2str((1-abs(MR_Signal_Pixel))*100) '%']);


%% Simulation functions
function [grad]=generate_gradient_points_local( seq , World , sys)

% Define first lobe
N_plat = ceil( seq.delta / sys.dt );
N_delay = ceil( ( seq.Delta - seq.delta ) / sys.dt );
Gx_pos = [0;  ones(N_plat , 1); 0];
Gx_delay = [0; zeros(N_delay , 1); 0];
Gx_neg = [0; -ones(N_plat , 1); 0];

grad = [Gx_pos; Gx_delay; Gx_neg ];

% grad = zeros( 2 + round( ( seq.delta+seq.Delta ) / sys.dt ) , 1 );
% grad( 1+1 : round( (seq.delta) / sys.dt) ) = 1;
% grad( 1 + round( seq.Delta / sys.dt ) : end-1 ) = -1;

% For now we have a normalized gradient shape
% We are going to calculate what is correct G_max of this
% waveform in mT/m to match the given b-value
Bfact=( (2 * pi * sys.gamma_bar ) ^ 2 ) * ( seq.delta ^ 2 ) *( seq.Delta - seq.delta / 3 );
G_max = sqrt( seq.bvalue / Bfact );

disp(['Target b-value = ',num2str(seq.bvalue./1e6),' s/mm^2']);
disp(['Design b-value = ' num2str( ( (2 * pi * sys.gamma_bar ) ^ 2 ) * (seq.delta^2) * (G_max^2) *(seq.Delta-seq.delta/3)*1e-6) ' s/mm^2' ] );
grad = grad .* G_max * 1e3;  % Gradient in mT/m
end

function [Dstep] = diffusion_step_local(diff)
Dstep=zeros(diff.N_spins,diff.dim);
for cpt=1:diff.dim
    KD = sqrt(diff.D(cpt) * 2 * diff.dim * diff.dt);
    Dstep(:,cpt)=KD * randn(diff.N_spins,1);
end
end

%% Display functions
function figH=create_figure_local(MR_Phase_Molecule,spin,Grad,cpt,Pixel,diff)

HandleFig = figure('Units', 'normal', 'Position', [0.35 0.1 0.55 0.8]); %[left bottom width height]

subgroup1 = axes(   'Parent', HandleFig, 'Units', 'normal', 'Position',  [0.0 0.3 0.4 0.4]);  % spin Displacement
subgroup2 = uipanel('Parent', HandleFig, 'Units', 'normal', 'Position',  [0.5 0.0 0.5 0.5]);  % Sequence
subgroup3 = axes(   'Parent', HandleFig, 'Units', 'normal', 'Position',  [0.55 0.6 0.45 0.4]);  % Phase diagram
% subgroup3 = axes(   'Parent', HandleFig, 'Units', 'normal', 'Position',  [0.55 0.1 0.45 0.8]);  % Phase diagram
subgroup4_X = axes(   'Parent', HandleFig, 'Units', 'normal', 'Position',[0.0 0.7 0.4 0.1]);  % Hist X
subgroup4_Y = axes( 'Parent', HandleFig, 'Units', 'normal', 'Position',  [0.4 0.3 0.1 0.4]);  % Hist Y


%% Spin Displacement
hold(subgroup1,'on');
% subgroup1.XLimMode='manual';
% subgroup1.YLimMode='manual';
figH.wp=plot(squeeze(spin(:,1,1:cpt))',squeeze(spin(:,2,1:cpt))','Parent',subgroup1);
figH.ws=scatter(squeeze(spin(:,1,cpt)),squeeze(spin(:,2,cpt)),50,'b','filled','Parent',subgroup1);

set(subgroup1,'Visible','off');
set(subgroup1,'xtick',[],'ytick',[]);
set(subgroup1,'Color','k');
set(gcf,'Color','k');

% subgroup1.XLimMode='manual';
% subgroup1.XLim=[-Pixel(1) Pixel(1)];
% subgroup1.YLimMode='manual';
% subgroup1.YLim=[-Pixel(2) Pixel(2)];


%% Sequence X
subgroup2.ForegroundColor='k';
subgroup2.BackgroundColor='k';
subgroup2.HighlightColor='k';
subgroup2.ShadowColor='k';
subgroup2_x = axes('Parent', subgroup2, 'Units', 'normal', 'Position', [0.1 0.55 0.9 0.40]);
hold(subgroup2_x,'on');
figH.g1=plot(1:cpt,zeros(cpt,1),'r','LineWidth',4,'Parent',subgroup2_x);

set(subgroup2_x,'box','off');
set(subgroup2_x,'Visible','on');
set(subgroup2_x,'xtick',[]);
% set(subgroup2_x,'ytick',[0 25 50 75 100]);
% set(subgroup2_x,'YTickLabel',[0 25 50 75 100]);
% set(subgroup2_x,'YMinorTick','on')
subgroup2_x.XAxis.Color = 'k';
subgroup2_x.YAxis.Color = 'w';
subgroup2_x.YAxisLocation = 'left';
set(subgroup2_x,'Color','w');
set(gcf,'Color','k');

subgroup2_x.XLimMode='manual';
subgroup2_x.XLim=[0-1000 size(Grad,1)+1000];
% subgroup2_x.YLimMode='manual';
% subgroup2_x.YLim=[0 120];
subgroup2_x.Color='k';
xtickangle(subgroup2_x,30);
xlabel(subgroup4_X,'X axis mT/m','Color','w')

%% Sequence Y
subgroup2_y = axes('Parent', subgroup2, 'Units', 'normal', 'Position', [0.1 0 0.9 0.40]);
hold(subgroup2_y,'on');
figH.g2=plot(1:cpt,zeros(cpt,1),'r','LineWidth',4,'Parent',subgroup2_y);
set(subgroup2_y,'Visible','on');
set(subgroup2_y,'xtick',[]);
% set(subgroup2_y,'ytick',[0 25 50 75 100]);
% set(subgroup2_y,'YTickLabel',[0 25 50 75 100]);
% set(subgroup2_y,'YMinorTick','on')
subgroup2_y.XAxis.Color = 'k';
subgroup2_y.YAxis.Color = 'w';
subgroup2_y.YAxisLocation = 'left';
set(subgroup2_y,'Color','w');
set(subgroup2_y,'box','off');

% set(gcf,'Color','k');

subgroup2_y.XLimMode='manual';
subgroup2_y.XLim=[0-1000 size(Grad,1)+1000];
% subgroup2_y.YLimMode='manual';
% subgroup2_y.YLim=[0 120];
subgroup2_y.Color='k';
xtickangle(subgroup2_y,15);
title(subgroup4_X,'Relative displacement','FontSize',12,'FontWeight','bold','Color','w')
xlabel(subgroup4_X,'X axis mT/m','Color','w')

%% Hist in X
RelativeDistanceX=squeeze(spin(:,1,cpt))-squeeze(spin(:,1,1));
figH.hx=histogram(RelativeDistanceX,20,'BinLimits',[-Pixel(1)/4 Pixel(1)/4],'Parent',subgroup4_X);

set(subgroup4_X,'ytick',[-Pixel(1)/5 0 Pixel(1)/5]);
set(subgroup4_X,'XTickLabel',[-1000*Pixel(1)/5 0 1000*Pixel(1)/5]);
set(subgroup4_X,'ytick',[]);
set(subgroup4_X,'box','off');
set(subgroup4_X,'Color','k');
set(subgroup4_X,'XMinorTick','on')
subgroup4_X.XAxis.Color = 'w';
subgroup4_X.YAxis.Color = 'k';
subgroup4_X.XLimMode='manual';
subgroup4_X.XLim=[-Pixel(1)/4 Pixel(1)/4];
subgroup4_X.YLimMode='manual';
subgroup4_X.YLim=[0 size(spin,1)*1/4];

title(subgroup4_X,'Relative displacement','FontSize',12,'FontWeight','bold','Color','w')
xlabel(subgroup4_X,'um','Color','w')

%% Hist in Y
RelativeDistanceY=squeeze(spin(:,2,cpt))-squeeze(spin(:,2,1));
figH.hy=histogram(RelativeDistanceY,20,'BinLimits',[-Pixel(2)/4 Pixel(2)/4],'Parent',subgroup4_Y);

set(subgroup4_Y,'ytick',[-Pixel(2)/5 0 Pixel(2)/5]);
set(subgroup4_Y,'XTickLabel',[-1000*Pixel(2)/5 0 1000*Pixel(2)/5]);
set(subgroup4_Y,'ytick',[]);
set(subgroup4_Y,'box','off');

set(subgroup4_Y,'Color','k');
set(subgroup4_Y,'View',[90 90]);
set(subgroup4_Y,'XMinorTick','on');
subgroup4_Y.XAxis.Color = 'w';
subgroup4_Y.YAxis.Color = 'k';

subgroup4_Y.XLimMode='manual';
subgroup4_Y.XLim=[-Pixel(2)/4 Pixel(2)/4];
subgroup4_Y.YLimMode='manual';
subgroup4_Y.YLim=[0 size(spin,1)*1/4];
xlabel(subgroup4_Y,'um','Color','w','Rotation',-90,'VerticalAlignment','top')
xtickangle(subgroup4_Y,-90)


set(gcf,'Color','k');

%% Spin Phase
hold(subgroup3,'on');
figH.pph=plot(1:cpt,MR_Phase_Molecule(1:end,1:cpt),'LineWidth',4,'Parent',subgroup3);

set(subgroup3,'Visible','on');
set(subgroup3,'Color','k');
%     set(subgroup3,'ytick',[-45 0 45]);
%     set(subgroup3,'yTickLabel',[-45 0 45]);
set(subgroup3,'xtick',[]);
set(subgroup3,'box','off');

set(subgroup3,'Color','k');
set(subgroup3,'YMinorTick','on');
subgroup3.XAxis.Color = 'k';
subgroup3.YAxis.Color = 'w';
subgroup3.XLimMode='manual';
subgroup3.XLim=[-1000 size(Grad,1)+1000];
subgroup3.YLimMode='manual';
subgroup3.YLim=[1.2*min(min(MR_Phase_Molecule)) 1.2*max(max(MR_Phase_Molecule))];
xlabel(subgroup3,'Phase (degree)','Color','w')
set(gcf,'Color','k');
end

function update_figure_local(figH,MR_Phase_Molecule,spin,Grad,cpt)

for cpt2=1:size(spin,1)
  figH.wp(cpt2).XData=squeeze(spin(cpt2,1,1:cpt))';
  figH.wp(cpt2).YData=squeeze(spin(cpt2,2,1:cpt))';
end

figH.ws.XData=squeeze(spin(:,1,cpt))';
figH.ws.YData=squeeze(spin(:,2,cpt))';

figH.g1.XData=1:cpt;
figH.g1.YData=squeeze(Grad(1:cpt,1));

figH.g2.XData=1:cpt;
figH.g2.YData=squeeze(Grad(1:cpt,2));

cpt_s=1;
for cpt2=1:size(spin,1)
  figH.pph(cpt_s).XData=1:cpt;
  figH.pph(cpt_s).YData=MR_Phase_Molecule(cpt2,1:cpt);
  cpt_s=cpt_s+1;
end

RelativeDistanceX=squeeze(spin(:,1,cpt))-squeeze(spin(:,1,1));
RelativeDistanceY=squeeze(spin(:,2,cpt))-squeeze(spin(:,2,1));
figH.hx.Data=RelativeDistanceX;
figH.hy.Data=RelativeDistanceY;

drawnow;
end