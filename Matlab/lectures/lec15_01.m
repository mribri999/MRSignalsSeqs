%
%	Model diffusion using EPG

T = 0.025; 			% seconds, gradient duration
A = 40;				% mT/m max gradient
kg = 2*pi*42.58*1000*A*T;
D = 1e-9*[0:0.1:3];			% diffusivity (m^2/s)
T1 = 3;
T2 = 2;					% seconds, make it long!
for d=1:length(D)
  Q = epg_m0;
  Q = epg_rf(Q,pi/2,pi/2);			% Excitation
  [Q,EE,B] = epg_grelax(Q,T1,T2,T,kg,D(d),1);	% Gradient
						% Note B is b-value for 1 grad.	
  Q = epg_rf(Q,pi,0);				% Refocusing pulse
  Q = epg_grelax(Q,T1,T2,T,kg,D(d),1);	% Gradient
  sig(d) = abs(Q(1,1));
end;

% B value from standard equation
b = (2*pi*42.58*1000*A*T)^2*T*(1-1/3);  
		% Assume gradient width and separation same!

plot(D,sig,'b--'); title('Signal vs Diffusivity');
hold on;
plot(D,exp(-b*D),'r:'); hold off;
legend('EPG Signal','Calculated');
xlabel('Diffusivity (m^2/s)');
ylabel('Signal');
grid on;



