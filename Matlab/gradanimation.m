% Example of an animation in Matlab

N = 24;			% number of frames
z = [-1:.05:1];		% number of locations

global framenum;	% declare as global
global filestem;	
framenum=0;		% Start at frame 0
filestem = '/Users/brian/tmp/ex'; 	% will save exNNNN.tif

for k=0:N
  M = [cos(pi*z*k/N); sin(pi*z*k/N); 0*z];	% twist vs n,z
  %plotm(M,1,[0*z;0*z;z]);			% could do 2x2
  showspins(M,1,[0*z;0*z;z]);			% one 3D plot
  saveframe;					% save the frame
end;

