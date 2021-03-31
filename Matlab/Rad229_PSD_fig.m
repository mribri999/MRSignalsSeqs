% This function generates a pulse sequence diagram plot.
%
% SYNTAX: ax = PAM_PSD_fig(RF, Gx, Gy, Gz, t)
%
% INPUTS: RF = B1 amplitude [1 x n] [microTesla]
%         Gx = Gx amplitude [t x n] [T/m]
%         Gy = Gy amplitude [t x n] [T/m]
%         Gz = Gz amplitude [t x n] [T/m]
%          t = time points  [t x 1] [s]
%
% OUTPUTS: ax  =  vector of graphics handles
%
% DBE@STANFORD.EDU (March 2021) for Rad229

function ax  =  PAM_PSD_fig(RF, Gx ,Gy, Gz, t)

if length(t)==1, t = ( 0 : length(RF)-1 ) * t; end % Expand t (i.e. dt) to a time vector
  
% Find the range of each input
bound = 0.02;
precision = 100;
RF_rng = max(RF)-min(RF); RF_min = min(RF)-0.1*RF_rng; RF_max = max(RF)+0.1*RF_rng;
Gx_rng = max(Gx(:))-min(Gx(:)); Gx_min = floor(precision*(min(Gx(:))-0.01*Gx_rng))/precision; Gx_max = ceil(precision*(max(Gx(:))+0.01*Gx_rng))/precision; 
  if Gx_min>-bound, Gx_min = -bound; end  % Define a maximum lower-bound 
  if Gx_max<+bound, Gx_max = +bound; end  % Define a minimum upper-bound
Gy_rng = max(Gy(:))-min(Gy(:)); Gy_min = floor(precision*(min(Gy(:))-0.01*Gy_rng))/precision; Gy_max = ceil(precision*(max(Gy(:))+0.01*Gy_rng))/precision;
  if Gy_min>-bound, Gy_min = -bound; end 
  if Gy_max<+bound, Gy_max = +bound; end
Gz_rng = max(Gz(:))-min(Gz(:)); Gz_min = floor(precision*(min(Gz(:))-0.01*Gz_rng))/precision; Gz_max = ceil(precision*(max(Gz(:))+0.01*Gz_rng))/precision;
  if Gz_min>-bound, Gz_min = -bound; end 
  if Gz_max<+bound, Gz_max = +bound; end

t_min = min(t)-0.02e-3;
t_max = max(t)+0.02e-3;

f = figure;
% RF axis
ax(1) = smart_subplot(4,1,1,0.05,0.05,0.1); hold on;
  pp(1) = plot([t_min t_max],[0 0],'k');
  p(1) = plot(t,RF);
  axis([t_min t_max RF_min-eps RF_max+eps]);
  ylabel('B_1 [µT]');

% X-gradient axis
ax(2) = smart_subplot(4,1,2,0.05,0.05,0.1); hold on;
  pp(2) = plot([t_min t_max],[0 0],'k');
  p(2) = plot(t,Gx);
  for n = 2:size(Gx,2)
    px(n-1) = plot(t,Gx(:,n));
  end
  axis([t_min t_max Gx_min-eps Gx_max+eps]);
  ylabel('G_x [T/m]');

% Y-gradient axis
ax(3) = smart_subplot(4,1,3,0.05,0.05,0.1); hold on;
  pp(3) = plot([t_min t_max],[0 0],'k');
  p(3) = plot(t,Gy(:,1));
  for n = 2:size(Gy,2)
    py(n-1) = plot(t,Gy(:,n));
  end
  axis([t_min t_max Gy_min-eps Gy_max+eps]);
  ylabel('G_y [T/m]');

% Z-gradient axis
ax(4) = smart_subplot(4,1,4,0.05,0.05,0.1); hold on;
  pp(4) = plot([t_min t_max],[0 0],'k');
  p(4) = plot(t,Gz);
  for n = 2:size(Gz,2)
    pz(n-1) = plot(t,Gz(:,n));
  end
  axis([t_min t_max Gz_min-eps Gz_max+eps]);
  ylabel('G_z [T/m]');
  xlabel('Time [s]');
    
% Stylize the figure
Rad229_plot_style(p);
% Rad229_fig_style(f);

% Set the colors for the zero-lines and phase encode lines
set(pp,'Color',[0.7 0.7 0.7]);

if size(Gx,2)>1
  set(px,'Color',[0.7 0.7 0.7]);
  set(px(end),'Color',get(p(2),'color'),'LineWidth',2);
end

if size(Gy,2)>1
  set(py,'Color',[0.7 0.7 0.7]);
  set(py(end),'Color',get(p(3),'color'),'LineWidth',2);
end

if size(Gz,2)>1
  set(pz,'Color',[0.7 0.7 0.7]);
  set(pz(end),'Color',get(p(4),'color'),'LineWidth',2);
end

return