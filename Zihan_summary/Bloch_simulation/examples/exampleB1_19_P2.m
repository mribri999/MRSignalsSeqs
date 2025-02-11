% Example with real RF pulse and steady state.
% Same RF, now excitation/recovery sequence.
%
% Starts the same as prior example...
% trying to calculate Ernst angle
clear;

%% RF pulse simulation
dt = .004;		% ms, sample spacing
% tip = 30;		% desired tip angle
tip = 1:30;		% desired tip angle (degree)
t = [-5:dt:5];		% extended time period to get 100 Hz spectral res.
tplot = find(abs(t)<=1); % only work on central 2 ms	
% rf=0*t;
rf=zeros(length(t),length(tip)); % allocate RF (related to t & different tip angle)
rfp = zeros(size(rf));
f = zeros(size(rf));
for i = 1:size(rf,2)
    rf(tplot,i) = msinc(length(tplot),3);	
    rf(:,i) = rf(:,i)/(sum(rf(:,i))*dt)*(tip(i)/360)/42.58;	% scale for desired flip.
    
    rfp(:,i) = fftshift(fft(fftshift(rf(:,i))))*dt*42.58;		% RF profile, scaled
    f(:,i) = ([1:length(rfp(:,i))]-length(rfp(:,i))/2)/length(rfp(:,i))/dt;	% kHz
end 

figure;		% Plot the RF and profile.
alpha2show = 10;
idx_alpha2show = find(tip==alpha2show);

subplot(2,1,1);
plot(t(tplot),rf(tplot,idx_alpha2show)); xlabel('time(ms)'); ylabel('B1 (mT)');title(['B1(t) for alpha = ',num2str(alpha2show),' (deg)']);
subplot(2,1,2);
freqplot = find(abs(f(:,idx_alpha2show))<5);			% Plot only range of freqs.
plot(f(freqplot,idx_alpha2show),abs(rfp(freqplot,idx_alpha2show))*360);	
title('Small Tip Approximation'); ylabel('Flip (deg)'); xlabel('Freq (kHz)');title(['slice-selective RF pulse profile (in freq) for alpha =',num2str(alpha2show),' (deg)'])

%% Sequence simulation (Hard pulse simulation)

TR = 5;			% ms
T1 = 500;		% ms
T2 = 50;		% ms

% Bloch simulation

Gz = 2.3;	% mT/m  (gamma/2pi*Gz ~ 100 kHz/m, 1kHz/cm; (Gz defined for slice-selective RF pulse)
df = 100e-3;		% kHz, off-resonance.

pos = [-.05:.0001:.05];	% 	Positions to simulate (for slice-selective RF pulse)

M = ones(3,length(tplot),length(pos),length(tip)); % simulate for all tip angle
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
for angle = 1:length(tip)
    for z = 1:length(pos)
      % tt = sprintf('Position %d of %d',z,length(pos)); disp(tt); % Show progress
      % Gradient rotation same for each interval
      Rgrad = zrot((42.58*pos(z)*Gz+df)*dt*360);	% zrot(a) = rotation matrix.
      A = As;	% Set to start value again! (relaxation over TR-RF duration & after perfect spoiling)
      B = Bs;	% Set to start value again. (relaxation over TR-RF duration & after perfect spoiling)

      for ti = 2:length(tplot)
        % Hard Pulse Approximation...
        alpha = rf(tplot(ti),angle)*dt*42.58*360;		% RF rotation over interval (simulate for that tip)
        Rrf = yrot(alpha);

        M(:,ti,z,angle) = Rrf*Rgrad*M(:,ti-1,z,angle);		% Apply RF, Gradient

        Ai = At*Rgrad*Rrf; % taking account in relaxation during RF
        Bi = Bt;

        % Propagate A,B
        A = Ai*A; % propagate A during RF pulse (accumulate based on As & Bs)
        B = Ai*B+Bi; % propagate A during RF pulse (accumulate based on As & Bs)

      end
      % Here A & B is the overall A & B for a whole TR 
      % - relaxation during TR-RF duration
      % - perfect spoiling after TR-RF duration
      % - zrot during slice selective gradient (i.e., Gz) (hard pulse approximation)
      % - yrot during RF pulse (hard pulse approximation)
      % - relaxation during RF pulse (hard pulse approximation)

      MM(:,z,angle) = inv(eye(3)-A)*B;	% Steady-state magnetization (matrix propagation)
      MMe(:,z,angle) = A*[0;0;1]+B;	% Magnetization after RF.  
                    % Note that
                    % if we start with M=M0 and propagate
                    % through 1 TR, there's no recovery before
                    % the RF, so we just get the excitation profile.

    end
end

%% Calc Ernst angle
centre_pos = floor(length(pos)/2)+1; % centre of pos
Mxy_centre = squeeze(MM(1,centre_pos,:) + i*MM(2,centre_pos,:)); % Mxy of steady-state at centre_pos with differnt tip angle
[Mxy_max,idx_ernst] = max(Mxy_centre);
disp(['Ernst angle is ',num2str(tip(idx_ernst)),' (deg)'])

figure;
plot(tip,Mxy_centre);
xlabel('Flip angle (deg)'), ylabel('Mxy of steady-state')
hold on
plot(tip(idx_ernst),Mxy_max,'ro','MarkerSize',8)
text(tip(idx_ernst),Mxy_max,sprintf('Ernst angle = %d (deg)',tip(idx_ernst)),'VerticalAlignment','top')
title('Flip angle & related steady state signal')


