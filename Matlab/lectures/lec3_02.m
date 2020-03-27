
% Lecture 3, Example 02
% 
% Use the rotation functions to animate excitation, precession
%
M = [0;0;1];	% Equilibrium
Ne = 10;		% Frames - excitation
Np = 50;		% Frames - precession
Nrot = 1;		% cycles
Ry = yrot(-90/Ne);	% Rotation about My
Rz = zrot(360/Nrot/Np);	% Rotation about Mz

figure(3); disp('Excitation/Rotation');
plotm(M);
for n=1:Ne; 
  M=Ry*M; plotm(M);	% Rotate and plot
  drawnow;
end;

for n=1:Np; 
  M=Rz*M; plotm(M);	% Rotate and plot
  drawnow;
end;

% Exercise:  Add relaxation so that the spins return to Mz after 5 rotations
