function h= imagescn(I, scale, dims, FigureWidth, border)
% function IMAGESCN(I,[min max],[rows cols],FigureWidth)
%
% function to display multiple images
%	I can be 2-d, i.e., single image
%            3-d,       array of images
%            4-d,       2-d array of images
%
% user specified scale [min max] is applied to all images
% NO--defaults to independent scaling of each sub-image to full range (for scale=[])
%     defaults to using range (min/max) of the entire input array (I) (DBE)
%
% sub-images are displayed in N1 rows by N2 columns
%   N1 and N2 may be specified in function input as "rows" and "cols"
%		if input values for rows & cols input are such that
%		rows*cols < full number of images, then the 1-st rows*cols images are displayed
%	or N1 and N2 may be calculated by default:
%		for 3-d array (N1 and N2 are calculated to make approx square image array)
%		for 4-d array (N1 and N2 are 3rd and 4th array dimensions, respectively)
%
% FigureWidth sets the width in inches (defaults to 6 inches). It also sets the paper width
% so that exported figures (to files) will have the specified width.
%
% usage:  imagescn(I), imagescn(I,scale), imagescn(I,[],[rows cols]), imagescn(I,scale,[rows cols])
%		  (note FigureWidth is optional parameter with any usage above)
%
% written by: 	Peter Kellman  (kellman@nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%				updated 07-Sep-2000 to fix bug
%

% Added by DBE
if nargin<=4
  border=0;
end

Nd=ndims(I);
if Nd==2; % case of single image
  N=1;
  N1=1; N2=1;
elseif Nd==3; % case of array of images
  N=size(I,3);
  N2=ceil(sqrt(N)); N1=ceil(N/N2);
elseif Nd==4; % case of 2-d array of images
%   N1=size(I,3);
%   N2=size(I,4);
  N1=size(I,4);
  N2=size(I,3);
  N=N1*N2;
end
if exist('dims');
  if length(dims)==2; rows=dims(1);cols=dims(2);
    N1=rows;N2=cols;
  else
    if ~isempty(dims);disp('Error: must enter [rows cols] for dimensions'); return;end
  end
end
% if ~exist('scale'); scale=[];end
if ~exist('scale')
  if min(I(:))~=max(I(:))
    scale=[min(I(:))-eps max(I(:))+eps];
  else
    scale=[min(I(:))-0.1*min(I(:))-eps max(I(:))+0.1*min(I(:))+eps];
  end
end  % Changed by DBE...

set(0,'Units','Inches');
scnsize=get(0,'ScreenSize'); % [left,bottom,width,height] % full screen size available
ScreenWidth=scnsize(3); ScreenHeight=scnsize(4); % width & height in inches
Xsize=size(I,2);Ysize=size(I,1); % size of pixels in image (Xsize x Ysize)
border_percent=.03;
deltaX =(border_percent*Xsize); deltaY=(border_percent*Ysize); % calculate pixel borders as % of size
X0size=N2*Xsize+(N2-1)*deltaX; Y0size=N1*Ysize+(N1-1)*deltaY; % full figure size in pixels (before display)
aspect_ratio=Y0size/X0size; % aspect ratio

% center figure on screen with specified figure width in inches
if ~exist('FigureWidth'); FigureWidth=6;end  % figure width in inches (default is 6")
FigureHeight=FigureWidth*aspect_ratio; % figure height in inches

% Added by DBE...ensures that entire figure is shown on screen;
if FigureHeight>ScreenHeight
  FigureWidth=0.9*FigureWidth*(ScreenHeight/FigureHeight);
  FigureHeight=0.9*ScreenHeight;
elseif FigureWidth>ScreenWidth
  FigureHeight=0.9*FigureHeight*(ScreenWidth/FigureWidth);
  FigureWidth=0.9*ScreenWidth;
end

FigureBottom=(ScreenHeight-FigureHeight)/2;
FigureLeft=(ScreenWidth-FigureWidth)/2;
set(gcf,'Units','Inches')
set(gcf,'Position',[FigureLeft FigureBottom FigureWidth FigureHeight])

% calculate sub-image dimensions in inches
SubImageWidth=FigureWidth*Xsize/X0size;
SubImageHeight=FigureHeight*Ysize/Y0size;
Xborder=FigureWidth*deltaX/X0size;
Yborder=FigureHeight*deltaY/Y0size;

% set background color to be white
set(gcf,'Color',[1 1 1]);

% calculate sub-image dimensions in normalized units
SubImageWidth=SubImageWidth/FigureWidth;
SubImageHeight=SubImageHeight/FigureHeight;
Xborder=Xborder/FigureWidth;
Yborder=Yborder/FigureHeight;

for k=1:min(N,N1*N2)
  i=ceil(k/N2);
  j=k-(i-1)*N2;
  SubImageLeft=(j-1)*(SubImageWidth+Xborder);
  SubImageBottom=(N1-i)*(SubImageHeight+Yborder);
  % 	subplot('Position',[SubImageLeft SubImageBottom SubImageWidth SubImageHeight])
  ax(i,j)=axes('Position',[SubImageLeft SubImageBottom SubImageWidth SubImageHeight]); % DBE
  if isempty(scale)
    index(i,j)  = imagesc(I(:,:,k)); axis image; axis off;
  else
    index(i,j)  = imagesc(I(:,:,k),scale); axis image;axis off
  end

  if border
    offset=1;
    x_pts=[0+offset size(I(:,:,k),2)-offset size(I(:,:,k),2)-offset 0+offset 0+offset];
    y_pts=[0+offset 0+offset size(I(:,:,k),1)-offset size(I(:,:,k),1)-offset 0+offset];
    l=line(x_pts,y_pts);
    set(l,'color',[0 0 0],'linewidth',1);
  end

end

set(gcf, 'PaperPosition', [0 0 FigureWidth FigureHeight]);

if nargout > 0
  %     h = index;
  h = ax;  % DBE
end;

% launch window level tool (button)
%wl_tool;
%pz_tool;
%roi_tool;