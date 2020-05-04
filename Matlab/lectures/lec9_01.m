
% Lecture 9, Example 01
% 
% Balanced SSFP signal plot

TR = 5;		% ms
TE = 0;		% ms
flip = 60;	% degrees
T1 = 1000;	% ms
T2 = 200;	% ms
rfdphase=180;    % degrees per TR

tt=sprintf('TR/TE = %d/%d ms, Flip=%d deg, T1/T2=%d/%d ms',TR,TE,flip,T1,T2);
disp(tt);

df = -200:200;
for k=1:length(df);
  s(k) = bssfp(flip,TR,TE,T1,T2,df(k)+1000/TR*(rfdphase/360));
end;


set(0,'defaultAxesFontSize',14);	% Default font sizes
set(0, 'DefaultLineLineWidth', 2);	% Default line width
magphase(df,s);
subplot(2,1,1); title('bSSFP signal vs Frequency (Hz)');
subplot(2,1,2); xlabel('Frequency (Hz)');


