function epg_saveframe
%function epg_saveframe
% 
%	If global variables framenum and filestem are global and defined,
%	captures the plot axis and saves to a tiff file 'path/fileNNNN.tif'
%	so that	a movie can be made, say, with Quicktime Pro 7.0
%
%	To use this first initialize by typing
%		global framenum filestem
%		framenum = 0;
%		filestem = 'path/file'
%	

global framenum filestem	% Declare global to this function
  
    % -- If framenum and filestem are defined,
if (length('framenum')>=1 && length('filestem')>=1)

        f = getframe(gcf);
        [fim,map] = frame2im(f);
        fn = sprintf('%s_%04d.tif',filestem,framenum);
        imwrite(fim,fn,'TIFF');
        framenum = framenum+1;
end;




