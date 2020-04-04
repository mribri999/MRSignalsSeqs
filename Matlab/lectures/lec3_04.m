% Lecture 3 - Example 4
% 
% Simulate/Plot a gradient of twist:

% -- Initialize
z = [-1:0.05:1];                % Z locations
M = [ones(size(z)); 0*z; 0*z];  % Excited magnetization.

% -- Phase Twist
ph = exp(pi*i*z);               % One cycle of phase.
phmult = [ph; conj(ph); 0*ph];  % 3xN phase Multiplier

% -- Apply and Plot
M = mc2mr(phmult.*mr2mc(M));    % Convert to Mx+i*My, add phase, convert back
disp('One cycle of Twist');
plotm(M,1,[0*z;0*z;z]);         % Plot twist!



