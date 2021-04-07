% Lecture 3 - Example 7
%
% Simulate/Plot an RF excitation

% -- Initialize
z = [-1:0.05:1];                % Z locations
M = [0*z; 0*z; ones(size(z))];  % Equilibrium magnetization.
anim = 1;			% Set to 1 to "animate"
if (anim==1)
  global framenum;		% Setup to save frames
  global filestem;		
  framenum=0;
  filestem='/Users/brian/tmp/movie/im' 	% Make image sequence
end;

Nrf = 100;			% #points
TB = 4;				% Time x Bandwidth product	
rf = msinc(Nrf,TB/2);
rf = rf*90/sum(rf);		% RF sums to 90 degrees.
% -- Phase Twist per point
ph = exp(pi*i*z*TB*3/Nrf);	% Time x Bandwidth cycles over slice,
				% during RF duration, and show 3x slice
phmult = [ph; conj(ph); ones(size(ph))];  % 3xN phase Multiplier

T = [1 i 0;1 -i 0;0 0 1];	% Real to Complex M 
Txy = [1 0 0 ;0 1 0;0 0 0];	% extract Mxy for plots
Mc = mr2mc(M);			% complex M vector

% -- Simulate RF excitation
for k=1:Nrf
  R = T*yrot(-rf(k))*inv(T);	% Rotation from RF
  Mc = phmult.*(R*Mc);		% Apply RF and gradient rotations.
  Mxy = Txy*real(mc2mr(Mc));	% Extract Mxy for plot
  if (anim==1)
    plotm(Mxy,1,[0*z;0*z;z]);	% Plot it.
    drawnow;
    epg_saveframe;
  end;
end;

% -- Simulate gradient refocusing
for k=1:floor(Nrf/2)+4		% The "+4" is empirically determined!!
  Mc = conj(phmult).*Mc;	% Refocusing gradient rotation.
  Mxy = Txy*real(mc2mr(Mc));	% Extract Mxy for plot
  if (anim==1)
    plotm(Mxy,1,[0*z;0*z;z]);	% Plot it.
    drawnow;
    epg_saveframe;
  end;
end;
   
% -- Final plot    
plotm(Mxy,1,[0*z;0*z;z]);	% Plot it.
drawnow;







  
  
