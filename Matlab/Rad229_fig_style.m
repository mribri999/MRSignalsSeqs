% This function makes some stylistic figure adjustments.
%
% SYNTAX: PAM_fig_style(handle,aspect,style,invert)
%
% EXAMPLE: PAM_fig_style(gcf,'rect',[],0);
%
% DBE@STANFORD.EDU (March 2021) for Rad229

function Rad229_fig_style(handle,aspect,style,invert)

%% Define the defaults...
if nargin<1, handle=gcf;    end
if nargin<2, aspect='rect'; end
if nargin<3, style ='PSD';  end % Default is for a pulse sequence diagram (RF, Gx, Gy, Gz)
if nargin<4, invert=0;      end
switch aspect
  case 'rect'
    width=6.5;
    height=4;
  case 'square'
    width=4;
    height=4;
  case 'tall'
    width=4;
    height=6.5;
  case 'wide'
    width=6.5;
    height=2;
  otherwise
    width=6.5;
    height=4;
end

if invert==0
  bg_color =[0.9 0.9 0.9];
  txt_color=[0.0 0.0 0.0];
else
  bg_color =[0.0 0.0 0.0];
  txt_color=[0.9 0.9 0.9];  
end

font_name='Helvetica';
font_size=10;
% format{1}='%6.1f';
% format{2}='%6.1f';
% format{3}='%6.1f';
format{1}='%6.1e';
format{2}='%6.1f';
format{3}='%6.1e';

%% Define the FIGURE properties
set(handle,'Units','inches');
pos=get(handle,'position');

set(handle,'Position',[pos(1) pos(2) width height],...
           'PaperSize',[width height],'PaperPosition',[0 0 width height],...
           'InvertHardCopy','Off','MenuBar','None','Color',bg_color);
cameratoolbar('hide');

txt=findobj(handle,'Type','text');
for n=1:numel(txt)
%   set(txt(n),'Color',[1 1 1],'FontSize',font_size);
  set(txt(n),'Color',txt_color,'FontSize',font_size);
end

%% Define the AXES properties
a=findobj(handle,'Type','Axes');

txt=[];
lab=[];
ind=1;
for n=1:length(a)
  % Adjust the font for each axis
%   set(a(n),'FontName',font_name,'FontSize',font_size,'FontWeight','Bold','Color',bg_color);
  set(a(n),'FontName',font_name,'FontSize',font_size,'Color',bg_color);
  set(a(n),'XColor',txt_color, 'YColor',txt_color,'ZColor',txt_color);
  
  % Collect the handles for each TEXT or LABEL object within the axis
  txt=[txt; findobj(a(n),'type','text')];
  lab(ind)=get(a(n),'XLabel'); ind=ind+1;
  lab(ind)=get(a(n),'YLabel'); ind=ind+1;
  lab(ind)=get(a(n),'ZLabel'); ind=ind+1;
  lab(ind)=get(a(n),'Title');  ind=ind+1;

  if isfield(get(a(n)),'XTickLabel'), set(a(n),'XTickLabel',num2str(get(a(n),'XTick')',format{1}),'FontSize',font_size); end
  if isfield(get(a(n)),'YTickLabel'), set(a(n),'YTickLabel',num2str(get(a(n),'YTick')',format{2}),'FontSize',font_size); end
  if isfield(get(a(n)),'ZTickLabel'), set(a(n),'ZTickLabel',num2str(get(a(n),'ZTick')',format{3}),'FontSize',font_size); end

  if strcmp(style,'PSD')
    if n>1, set(a(n),'XTickLabel',[]); end % Turn off xlabels for RF, Gx, and Gy
  end
end

%% Set the defaults for the Labels and Titles
for n=1:length(txt)
%   set(txt(n),'FontName',font_name,'FontSize',font_size,'FontWeight','Bold','Color',txt_color);
  set(txt(n),'FontName',font_name,'FontSize',font_size,'Color',txt_color);
end

for n=1:length(lab)
%   set(lab(n),'FontName',font_name,'FontSize',1.25*font_size,'FontWeight','Bold','Color',txt_color);
  set(lab(n),'FontName',font_name,'FontSize',1.25*font_size,'Color',txt_color);
end

% % If there are a lot of labels, then rotate them for legibility
% if length(get(a(1),'XTickLabel'))>5
%   set(gca,'XTickLabelRotation',45);
% end
return