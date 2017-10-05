% Effect of Crushers

G = 40;		% mT/m;
T = 1.174;	% ms.
gamma = 42.58	% kHz/mT
alpha = [180:-1:60];

% 2 cycles/mm.


x = [0:.01:1]/1000;	% m.

Grot = 360*x*G*T*gamma; 	% Rotation due to gradient at each voxel.

Ms=[1;0;0];		% Mstart

for n=1:length(alpha)
 for k=1:length(x)
  M = Ms;
  M = zrot(Grot(k))*M;		% Crusher 1
  M = xrot(alpha(n))*M;		% Refocusing pulse.
  M = zrot(Grot(k))*M;		% Crusher 2
  Mend(:,k)=M;
 end;

 %figure(3);
 %plot(abs(Mend(1,:)+i*Mend(2,:)));
 figure(1);
 Mxy = Mend(1,:)+i*Mend(2,:);
 plot(x*1000,real(Mxy),'k--'); hold on;
 axis([min(x)*1000 max(x)*1000 -1.2 1.2]);
 plot([min(x) max(x)]*1000,[1 1]*abs(mean(Mxy)),'b-'); hold off;
 grid on;
 xlabel('Position (mm)'); ylabel('Signal');
 legend('M_{xy}','Avg M_{xy}');
 tt = sprintf('%d Degree Refocusing Angle',alpha(n)); title(tt);
 setprops;
 drawnow;
 %fig2tiff('crush',n);
 n
 Mse(n) = abs(mean(Mxy));
end;
 

figure(2);
plot(alpha,Mse); grid on; xlabel('Refocusing Angle (deg)');
ylabel('Spin Echo Signal'); title('Spin Echo vs Refoc. Angle');
a = axis; axis([a(1:2) 0 1]);
setprops;


