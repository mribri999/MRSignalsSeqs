% This function simply places the subplot axes with 
% better positioning than the default call to SUBPLOT.
% The axis handle is returned as 'h'.
%
% SYNTAX : h=smart_subplot(m,n,p[,gap,lm,bm])
%
% INPUTS : m   - Number of rows
%          n   - Number of columns
%          p   - Plot number (raster order)
%          gap - Gap between plots in normalized units [DEFAULT=0.05]
%          lm  - Left margin in normalized units [DEFAULT=0.0]
%          bm  - Bottom margin in normalized units [DEFAULT=0.0]
%
% OUTPUTS: h   - vector of handles for each of the subplot axes
%
% DBE 09/18/03

function h=smart_subplot(m,n,p,gap,left_margin,bot_margin)

if nargin<3
  error('Insufficient number of input arguments.');
end

if p>m*n
  error('Requested subplot # (p) was beyond range of m*n.');
end

% c_gap=0.05;
c_gap=0.00;

if nargin<4 | isempty(gap),         gap=0.05;      end
if nargin<5 | isempty(left_margin), left_margin=0; end
if nargin<6 | isempty(bot_margin),  bot_margin=0;  end

% if nargin==3
%   gap=0.1/min([m n]);
%     if gap>0.1, gap=0.1; end
% end

% The indices are used such that the produce the *same* order as SUBPLOT
[j,i]=ind2sub([n m],p);

width =((1-left_margin)-c_gap-gap*(n+1))/n;
height=((1-bot_margin)-gap*(m+1))/m;
left  =(j-(1-c_gap))*(gap+width)+gap+left_margin;
bottom=1-(i*(gap+height));

% [left bottom width height]
h=axes('position',[left bottom width height]);
% set(gca,'color',[1 1 1],'FontWeight','Bold','XColor',[1 1 1],'YColor',[1 1 1],'FontSize',15);

return

% function h=smart_subplot(m,n,p,gap)
% 
% if nargin<3
%   error('Insufficient number of input arguments.');
% end
% 
% if p>m*n
%   error('Requested subplot # (p) was beyond range of m*n.');
% end
% 
% % if nargin==3
% %   gap=0.1/min([m n]);
% %     if gap>0.1, gap=0.1; end
% % end
% 
% % The indices are used such that the produce the *same* order as SUBPLOT
% [j,i]=ind2sub([n m],p);
% 
% width =(1-gap*(n+1))/n;
% height=(1-gap*(m+1))/m;
% left  =(j-1)*(gap+width)+gap;
% bottom=1-(i*(gap+height));
% 
% h=subplot('position',[left bottom width height]);
% set(gca,'color',[1 1 1],'FontWeight','Bold','XColor',[1 1 1],'YColor',[1 1 1],'FontSize',15);
% 
% return