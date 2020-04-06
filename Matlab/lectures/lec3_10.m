% Lecture 3 - Example 10
%
% Simulate RF with gradient twist

TR = 0.01;		% s
T1 = 1;			% s
T2 = 0.2;		% s
% -- Initialize
z = [-1:0.05:1];                % Z locations
Mc = [0*z; 0*z; ones(size(z))];  % Equilibrium
[A,B] = relax(TR,T1,T2);	% Relaxation over TR
B = B*ones(size(z));		% Expand B vector 

Rex = yrot(-30);		% 30 degree excitation
T = [1 i 0; 1 -i 0; 0 0 1];	% Transform real to complex
Rex = T*Rex*inv(T);		% R in complex-m representation

% -- Phase Twist
ph = exp(pi*i*z);		% One cycle of phase.
phmult = [ph; conj(ph); ones(size(ph))];  % 3xN phase Multiplier

% -- Simulate to Steady State
figure(1);

for k=1:200			% Steady state so just start after RF
  Mc = A*Mc+B;			% Relaxation				
  Mc = phmult.*Mc;		% Gradient twist
  Mc = Rex * Mc;		% Excitation
  %plotm(real(mc2mr(Mc))); 
  %drawnow; pause;
end;

plotm(mc2mr(Mc),0.5);
 







  
  
