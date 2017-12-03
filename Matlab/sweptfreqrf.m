%
%%	Demo of swept frequency RF.
%
%	


zlocs = [-5:.01:5];		% cm
gamma = 4258;		% Hz/G
T = 3;			% ms;
a = 0.1;		% G

Gz = 3;			% G/cm
flim = gamma*Gz* [min(zlocs) max(zlocs)]/1000;	% Freq in kHz
freq = [-1:0.005:1] * flim(2);		% Sweep frequency linearly

t = [0:length(freq)-1]/length(freq)*T;			% ms
dT=t(2)-t(1);				% ms

b1 = a * ones(size(freq)) .* exp(2*pi*i*cumsum(freq)*dT);
gz = Gz * ones(size(freq));

b1 = [b1 0*b1];
gz = [gz -gz];

t = [t T+t];

figure(1);
magphase(t,b1);
subplot(2,1,1); title('RF pulse Magnitude');
setprops;


[mx,my,mz] = blochsim(b1,gz,dT,0,10*z);

M = [zeros(2,length(zlocs)); ones(1,length(zlocs))];

figure(2);

for n=1:length(b1)
  for z = 1:length(zlocs)
    M(:,z) = zrot((gamma*gz(n)*zlocs(z)/1000)*dT*360)* M(:,z); % Gradient
    M(:,z) = throt(gamma*abs(b1(n))/1000*dT*360,180/pi*angle(b1(n)))*M(:,z);  % RF
    mxy = M(1,:)+i*M(2,:);	% Transverse Magnetization
  end;

  vox =   ones(1,20)/20;
  mxys = conv2(mxy,vox,'same');	% Average over voxel.
 
  magphase(zlocs,M(1,:)+i*M(2,:));
  subplot(2,1,1); title('Transverse Magnetization');
  axis([min(zlocs) max(zlocs) 0 .5]);
  hold on;
  plot(zlocs,abs(mxys),'r--'); hold off;
  legend('M_{xy}','Voxel Signal');
  grid on;
  drawnow; 
  setprops;
  if (1==1)
    fname = sprintf('/Users/brian/tmp/sf/im.%04d.tif',n);
    print(gcf,'-dtiff',fname);
  end;
end;


