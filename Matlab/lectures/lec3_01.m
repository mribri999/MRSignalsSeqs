
% Lecture 3, Example 01
% 
% Plot some mangetization vectors using plotm.m
%
% Simple 3x1 vectors - ; ends a row in Matlab:
M1 = [0; 0; 1];         % "Equilibrim" magnetization along Mz
M2 = [1; 0; 0];         % Magnetization long Mx
M3 = [0; 0.5; 0];       % Shorter magnetization along My

% plotm.m can plot a (3xN) array of spins all at once:
figure(1); disp('Example of 3 Spins');
plotm([M1 M2 M3]);      % Notice we concatenate them horizontally.
%
% Which color corresponds to which M vector?
%
% You can also plot some dephased spins:
%
figure(2); disp('Dephased Spins');
x = [0:15];
M = [cos(2*pi*x/15); sin(2*pi*x/15); 0*x];
plotm(M);
