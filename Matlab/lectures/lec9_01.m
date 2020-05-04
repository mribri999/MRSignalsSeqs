% Lecture 9, Example 01
%
% Balanced SSFP signal plot

TR = 5;         % ms
TE = TR/2; % ms
flip = 30;      % degrees
T1 = 1000;      % ms
T2 = 200;       % ms
rfdphase=180;    % degrees per TR
flips = [60,90,30,10,5,1];    % degrees
fliplist = [1,4];             % Just do one flip angle to start
labels = {'90','60','30','10','5','1'};

tt=sprintf('TR/TE = %d/%d ms, Flip=%d deg, T1/T2=%d/%d ms',TR,TE,flip,T1,T2);
disp(tt);

df = -200:200;
s = zeros(length(df),length(flips(fliplist)));

for f=1:length(flips(fliplist))
  for k=1:length(df);
    s(k,f) = bssfp(flips(fliplist(f)),TR,TE,T1,T2,df(k)+1000/TR*(rfdphase/360));
  end;
end;

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width
figure(1);
subplot(2,1,1);
plot(df,abs(s)); 
lplot('Frequency(Hz)','Signal Magnitude','bSSFP Signal vs Frequency (Hz)');
legend(labels{fliplist});

subplot(2,1,2);
plot(df,angle(s)/pi); 
lplot('Frequency(Hz)','Signal Phase/\pi','bSSFP Signal vs Frequency (Hz)');
legend(labels{fliplist});

figure(2);
plot(real(s),imag(s),'.');
lplot('M_x','M_y','Signal in Mx-My Plane');
axis equal;

