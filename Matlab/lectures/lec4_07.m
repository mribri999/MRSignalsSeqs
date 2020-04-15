% Lecture 4, Example 07
% 
% Apply simple exp(i*phi) and see effect
Q = [0 1;0 1;0 0.5];	% one cycle in F+,F-,Z

if (1==1) % Setup for animation/movie capture
  global filestem;
  global framenum;
  filestem='/Users/brian/tmp/ph';
  framenum=0;
end;

Nf=120;
for k=1:Nf
  epg_show(Q*exp(2*pi*i*k/Nf),[],[],39,1); drawnow;
end;

