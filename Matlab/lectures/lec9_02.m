% Lecture 9, Example 02
%
% GRE signal, vs bSSFP.

TR = 5;         % ms
TE = 0; 	% ms
flip = 30;      % degrees
T1 = 1000;      % ms
T2 = 200;       % ms
rfdphase=180;    % degrees per TR
flips = [60,90,30,10,5,1];    % degrees
%fliplist = [1];             % Just do one flip angle to start
fliplist=[1 3];
labels = {'60','90','30','10','5','1'};

tt=sprintf('TR/TE = %d/%d ms, Flip=%d deg, T1/T2=%d/%d ms',TR,TE,flip,T1,T2);
disp(tt);

df = -200:200;
s = zeros(length(df),length(flips(fliplist)));
gre=s;

for f=1:length(flips(fliplist))
  for k=1:length(df);
    s(k,f) = bssfp(flips(fliplist(f)),TR,TE,T1,T2,df(k)+1000/TR*(rfdphase/360));
    gre(k,f) = gresignal(T1,T2,TE,TR,flips(fliplist(f)));
  end;
end;

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width
figure(1);
subplot(2,1,1);
plot(df,abs(s),'--');  hold on;
plot(df,abs(gre),'-');  hold off;
lplot('Frequency(Hz)','Signal Magnitude','bSSFP and Gradient-Spoiled Signals');
legend({labels{fliplist},labels{fliplist}});

subplot(2,1,2);
plot(df,angle(s)/pi,'--'); hold on;
plot(df,angle(gre)/pi,'-'); hold off;
lplot('Frequency(Hz)','Signal Phase/\pi','bSSFP Signal vs Frequency (Hz)');
legend({labels{fliplist},labels{fliplist}});

figure(2);
plot(real(s),imag(s),'.'); hold on;
plot(real(gre),imag(gre),'+'); hold off;
lplot('M_x','M_y','Signal in Mx-My Plane');
axis equal;
legend(labels{fliplist});

