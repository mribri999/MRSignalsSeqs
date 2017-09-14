%
%  Plot signal and EPG states for varying 90-180 phase.

etl = 64;
ph = [0:10:180];
TR = [1:etl];
flip = 40;

for k=1:length(ph)
  rf = [90+flip/2 ones(1,etl-1)*flip]*exp(i*ph(k)*pi/180)*pi/180;
  [s,ps] = epg_cpmg(rf ,[],1000,200,10);
  if (k==1) s0 = s; end;
  subplot(1,2,2);
  imagesc(abs(ps)); colormap('jet'); colorbar;
  xlabel('TR'); ylabel('F/Z states');

  subplot(1,2,1);
  plot(TR,abs(s0),'b:',TR,real(s),'b--',TR,imag(s),'r--',TR,abs(s),'k-');
  legend('CPMG','M_x','M_y','Abs(signal)');
  grid on; 
  xlabel('TR'); ylabel('Signal');
  tt = sprintf('%d deg Deviation from CPMG Angle',ph(k)); title(tt);
  a = axis; axis([0 etl -1 1]);
  %setprops;
  drawnow;
  %fig2tiff('cpmg',k); 	% Save figure for movie.
end;

