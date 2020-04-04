% Lecture 3 question 2

global framenum;
global filestem;
filestem='/Users/brian/tmp/q2'
framenum = 0;

[A,B] = relax(0.05,2,1);		% Relaxation
AA = A*xrot(30);		% Do rotation first.

M = [0;0;1];

for k=1:120
  plotm(M);
  M = AA*M+B;
  drawnow;
  saveframe
end;


