% Example with real RF pulse and steady state.
% Same RF, now excitation/recovery sequence.
%
% Starts the same as prior example...
clear;

%% RF pulse simulation
dt = .004;		% ms, sample spacing
tip = 30;		% desired tip angle
t = [-5:dt:5];		% extended time period to get 100 Hz spectral res.
tplot = find(abs(t)<=1); % only work on central 2 ms
rf=0*t;			% allocate RF
rf(tplot) = msinc(length(tplot),3);	% put RF pulse in central part.
rf = rf/(sum(rf)*dt)*(tip/360)/42.58;	% scale for desired flip.

rfp = fftshift(fft(fftshift(rf)))*dt*42.58;		% RF profile, scaled
f = ([1:length(rfp)]-length(rfp)/2)/length(rfp)/dt;	% kHz

figure(1);		% Plot the RF and profile.
subplot(2,1,1);
plot(t(tplot),rf(tplot)); xlabel('time(ms)'); ylabel('B1 (mT)');title('B1(t)');
subplot(2,1,2);
freqplot = find(abs(f)<5);			% Plot only range of freqs.
plot(f(freqplot),abs(rfp(freqplot))*360);	
title('Small Tip Approximation'); ylabel('Flip (deg)'); xlabel('Freq (kHz)');title('slice-selective RF pulse profile (in freq)')

%% Sequence simulation

TR = 5;			% ms
T1 = 500;		% ms
T2 = 50;		% ms

% Bloch simulation

Gz = 2.3;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm; (Gz defined for slice-selective RF pulse)
df = 100e-3;		% kHz, off-resonance.

pos = [-.05:.0001:.05];	% 	Positions to simulate (for slice-selective RF pulse)

M = ones(3,length(tplot),length(pos));
M(1:2,:)=0;				% M=[1;0;0]; the initial M

% -- Start right after RF pulse - here's where steady state will be!
% -- First set A,B for relaxation over (TR-RF duration)

A = diag([exp(-(TR-2)/T2) exp(-(TR-2)/T2) exp(-(TR-2)/T1)]); % (TR-2), according to the above RF pulse simulation, the duration of RF pulse is 2 ms
Bs = [0;0;1-exp(-(TR-2)/T1)];
% [A,Bs] = relax(TR-2,T1,T2,0); % This is an alternative way for the above two lines
As = [0 0 0;0 0 0;0 0 1]*A;	% Perfect Spoiling


% -- Relaxation during RF pulse:  A/B for each interval.
% considering relaxation during RF pulse
At = diag([exp(-(dt)/T2) exp(-(dt)/T2) exp(-(dt)/T1)]);
Bt = [0;0;1-exp(-(dt)/T1)];
% % do not considering relaxation during RF pulse
% At = diag([exp(-(0)/T2) exp(-(0)/T2) exp(-(0)/T1)]); % diag(1,1,1)
% Bt = [0;0;1-exp(-(0)/T1)]; % [0 0 0]


% -- Hard pulse approximation, with relaxation.
for z = 1:length(pos)
  tt = sprintf('Position %d of %d',z,length(pos)); disp(tt); % Show progress
  % Gradient rotation same for each interval
  Rgrad = zrot((42.58*pos(z)*Gz+df)*dt*360);	% zrot(a) = rotation matrix.
  A = As;	% Set to start value again! (relaxation over TR-RF duration & after perfect spoiling)
  B = Bs;	% Set to start value again. (relaxation over TR-RF duration & after perfect spoiling)
   
  for ti = 2:length(tplot)
    % Hard Pulse Approximation...
    alpha = rf(tplot(ti))*dt*42.58*360;		% RF rotation over interval
    Rrf = yrot(alpha);

    M(:,ti,z) = Rrf*Rgrad*M(:,ti-1,z);		% Apply RF, Gradient

    Ai = At*Rgrad*Rrf; % taking account in relaxation during RF
    Bi = Bt;

    % Propagate A,B
    A = Ai*A; % propagate A during RF pulse (accumulate based on As & Bs)
    B = Ai*B+Bi; % propagate A during RF pulse (accumulate based on As & Bs)
   
  end;
  % Here A & B is the overall A & B for a whole TR 
  % - relaxation during TR-RF duration
  % - perfect spoiling after TR-RF duration
  % - zrot during slice selective gradient (i.e., Gz) (hard pulse approximation)
  % - yrot during RF pulse (hard pulse approximation)
  % - relaxation during RF pulse (hard pulse approximation)

  MM(:,z) = inv(eye(3)-A)*B;	% Steady-state magnetization (matrix propagation)
  MMe(:,z) = A*[0;0;1]+B;	% Magnetization after RF.  
                % Note that
				% if we start with M=M0 and propagate
				% through 1 TR, there's no recovery before
				% the RF, so we just get the excitation profile.

end;

%% Display

Mxy = MM(1,:)+i*MM(2,:); % Mxy of steady-state
Mxye = MMe(1,:)+i*MMe(2,:); % Mxy of just after RF pulse (without relaxation after RF pulse)

figure(2);
plot(pos*100,abs(Mxy),'k-',pos*100,real(Mxy),'b-',pos*100,imag(Mxy),'r:');
grid on;
xlabel('Position (cm)');
ylabel('M_{xy}');legend('abs(Mxy)','real(Mxy)','imag(Mxy)');
title('Excitation/Recovery Signal (M for steady-state) vs Position');
setprops;

figure(3);
subplot(2,1,1);
plot(pos*100,abs(Mxye),'k-',pos*100,real(Mxye),'b-',pos*100,imag(Mxye),'r:');
grid on;
xlabel('Position (cm)');
ylabel('Me_{xy}');legend('abs(Mxye)','real(Mxye)','imag(Mxye)');
title('Excitation/Recovery Signal (M after excitation) vs Position');
setprops;

subplot(2,1,2);
plot(pos*100,180/pi*asin(abs(Mxye)));		% Plot of *ANGLE*
grid on;
xlabel('Position (cm)');
ylabel('Flip Angle (deg)');
title('Excitation Profile vs Position (B1 map); flip angle less than 30 deg since relaxation during RF');
setprops;

% ZN: I don't get this 
% figure(4);
% fplot('sin(pi/180*x)*(1-exp(-5/500))/(1-exp(-5/500)*cos(pi/180*x))',[0,60]);
% lplot('Signal','Flip Angle (deg)','Signal vs Flip Angle');
% setprops;

