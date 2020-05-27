% Lecture 14, Example 03
%
% Show M curves for different T1s, with inversion-recovery.

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width

T1 = [100 300 500 1000 1500 2000];
m = ones(size(T1));			% Note we only care about Mz!
invtimes = [0.5, 1000];		% Times where inversions played
t = [0:5000];				% Time-course in ms

E1 = exp(-(t(2)-t(1))./T1);


for tt = 1:length(t)-1
  f = find(invtimes>=t(tt) & invtimes<t(tt+1));
  if (length(f)>0)
    txt=sprintf('Inversion at time %d',t(tt)); disp(txt);
    m(tt+1,:) = -m(tt,:).*E1+(1-E1); 	% Inversion
  else
    m(tt+1,:) = m(tt,:).*E1+(1-E1); 	% No Inversion
  end;
end;

plot(t,m);
legend('100ms','300ms','500ms','1000ms','1500ms','2000ms');
lplot('Time (ms)','Signal','Signal vs Time for Multiple Inversion');
hold on;
plot(t,zeros(size(t)),'k-');	% add line along 0
hold off;

  
