%% Simulate RF pulse over different positions and time

dt = .004;		% ms, sample spacing
tip = 90;		% desired tip.	
t = [-5:dt:5];		% extended time period to get 100 Hz spectral res.
tplot = find(abs(t)<=1); % only work on central 2 ms
rf=0*t;			% allocate RF
rf(tplot) = msinc(length(tplot),3);	% put RF pulse in central part.
rf = rf/(sum(rf)*dt)*(tip/360)/42.58;	% scale for desired flip.

rfp = fftshift(fft(fftshift(rf)))*dt*42.58;		% RF profile, scaled
f = ([1:length(rfp)]-length(rfp)/2)/length(rfp)/dt;	% kHz

figure(1);		% Plot the RF and profile.
subplot(3,1,1);
plot(t(tplot),rf(tplot)); xlabel('time(ms)'); ylabel('B1 (mT)'); title('B1(t)')
subplot(3,1,2);
freqplot = find(abs(f)<5);			% Plot only range of freqs.
plot(f(freqplot),abs(rfp(freqplot))*360);	
title('Small Tip Approximation'); ylabel('Flip (deg)'); xlabel('Freq (kHz)');


%% Bloch simulation

Gz = 2.5;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm;
df = 300e-3;		% kHz, off-resonance.

pos = [-.05:.0001:.05];	% 	Positions to simulate

M = ones(3,length(tplot),length(pos));
M(1:2,:)=0;				% M=[1;0;0];

%plot Gz
subplot(3,1,3);
plot(pos*10,pos*Gz);
xlabel('Position(cm)'); ylabel('Delta_B0 (mT)'); title('Changed B0 along z axis')


% Note that we neglect relaxation during the RF.

for z = 1:length(pos)
  % Gradient rotation same for each interval
  Rgrad = zrot((42.58*pos(z)*Gz+df)*dt*360);	% zrot(a) = rotation matrix.
   
  for ti = 2:length(tplot)
    % Hard Pulse Approximation...
    alpha = rf(tplot(ti))*dt*42.58*360;		% RF rotation over interval
    Rrf = yrot(alpha);

    M(:,ti,z) = Rrf*Rgrad*M(:,ti-1,z);		% Apply RF, Gradient

  end;
end;

%% To plot
figure(2);
subplot(3,2,1);
plot(t(tplot),squeeze(M(1,:,ceil(length(pos)/2)))); 
xlabel('Time (ms)'); ylabel('M_x(t)'); title('M_x(t) in the centre of pos')

subplot(3,2,3);
plot(t(tplot),squeeze(M(2,:,ceil(length(pos)/2))));
xlabel('Time (ms)'); ylabel('M_y(t)'); title('M_y(t) in the centre of pos')

subplot(3,2,5);
plot(t(tplot),squeeze(M(3,:,ceil(length(pos)/2))));
xlabel('Time (ms)'); ylabel('M_z(t)'); title('M_z(t) in the centre of pos')

subplot(3,2,2)
plot(pos*100,squeeze(M(1,end,:)),'b--'); hold on;
plot(pos*100,squeeze(M(2,end,:)),'r-'); hold off;
xlabel('Position (cm)'); ylabel('M_x(z) and M_y(z)'); 
title('M_x(z) and M_y(z) at the end of time (final state)')
legend('M_x(z)','M_y(z)')

subplot(3,2,4)
plot(pos*100,abs(squeeze(M(1,end,:)+i*M(2,end,:))));
xlabel('Position (cm)'); ylabel('M_{xy}(z)');
title('M_{xy}(z) at the end of time (final state)')

subplot(3,2,6)
plot(pos*100,squeeze(M(3,end,:)));
xlabel('Position (cm)'); ylabel('M_y(z)');
title('M_z(z) at the end of time (final state)')

