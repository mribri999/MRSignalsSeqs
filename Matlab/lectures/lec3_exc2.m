% Exercise to make a Bloch simulator (space and time)

% function M = blochsim229(b1,gr,tp,t1,t2,df,dp)
%
%	Function does a Bloch simulation over time/space
%
%       INPUT:
%               b1 = (Mx1) RF pulse in mT.  Can be complex.
%               gr = (Mx1,2,or 3) 1,2 or 3-dimensional gradient in mT/m.
%               tp = (Mx1) time duration of each b1 and gr point, in seconds.
%               t1 = T1 relaxation time in seconds.
%               t2 = T2 relaxation time in seconds.
%               df = (Nx1) Array of off-resonance frequencies (Hz)
%               dp = (Px1,2,or 3) Array of spatial positions (cm).
%                       Width should match width of gr.
%
%	OUTPUT:
%		M = 3xPxN array of outputs


Mout = zeros(3,length(b1),length(dp));	% Allocate output

% Form E1, E2 and relaxation matrix/vector

for t=1:length(b1)

% Make RF rotation matrix for point

% Make Gradient rotation for point (loop if doing as real-valued)

  % Apply df rotation
  % Apply RF rotation
  % Apply gradient rotation
  % Apply relaxation

end;



