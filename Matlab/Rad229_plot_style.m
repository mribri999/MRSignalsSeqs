% This function defines a color palette for plots and can also be used to
% fetch the default colors for other uses (e.g. colormap_style.m).
%
% color = Rad229_plot_style(p)
%
% DBE@STANFORD.EDU (March 2021) for Rad229

function color = Rad229_plot_style(p)

% Color Pallette - In a preferred order...
if nargin==0 || numel(p)==3  % Permute the colors differently under the assumption that three plots is for Mx, My, Mz
  color_sat(1,:)=[1.0 0.0 0.0];  % Red
  color_sat(2,:)=[0.0 1.0 0.0];  % Green
  color_sat(3,:)=[0.0 0.0 1.0];  % Blue
  color_sat(4,:)=[1.0 0.5 0.0];  % Orange
  color_sat(5,:)=[1.0 1.0 0.0];  % Yellow
  color_sat(6,:)=[0.5 0.0 1.0];  % Purple
  color_sat(7,:)=[0.0 1.0 1.0];  % Cyan
elseif numel(p)>3  % Permute the colors differently under the assumption that three plots is for B1, Gx, Gy, Gz
  color_sat(1,:)=[1.0 0.5 0.0];  % Orange
  color_sat(2,:)=[1.0 0.0 0.0];  % Red
  color_sat(3,:)=[0.0 1.0 0.0];  % Green
  color_sat(4,:)=[0.0 0.0 1.0];  % Blue
  color_sat(5,:)=[1.0 1.0 0.0];  % Yellow
  color_sat(6,:)=[0.5 0.0 1.0];  % Purple
  color_sat(7,:)=[0.0 1.0 1.0];  % Cyan
else
  error('Input argument (number of plot handles) was incorrect.');
end

HSV=rgb2hsv(color_sat);

HSV(:,2)=0.50; % Reduce the saturation
HSV(:,3)=0.80; % Reduce the brightness

color=hsv2rgb(HSV);

% Normalized for luminance
M=repmat(mag(color,2),[1 3]);
color=color./M;

% Set the colors, if a PLOT handle was passed
if nargin>0
  for n=1:numel(p)
%     if isfield(p(n),'LineWidth'), set(p(n),'LineWidth',3);       end
%     if isfield(p(n),'Color'),     set(p(n),'Color',color(n,:));  end
    set(p(n),'LineWidth',3);
    set(p(n),'Color',color(n,:));
  end
end

return