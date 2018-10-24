%
%	Example of EPG animations (CPMG sequences)
%
%	Some parameters:  (Remove % as needed!)

%flips = [180 180 180 180 180] * pi/180;		% 180-refocusing angles
%flips = [150 120 120 120 120] * pi/180;		% 120-refocusing angles
%flips = [120 120 120 120 120] * pi/180;		% 120-refocusing angles
%flips = [150 120 120 120 120]*i * pi/180;	% 120-refoc. Non CPMG
flips = [120 60 60 60 60] * pi/180;		% 60-refocusing angles

% -- These lines setup a short hyperecho sequence:
%flips = [1 1.3 1.2 1.3] .* exp(pi*i*[0.62 -.28 .28 .25]);
%flips=[flips pi -fliplr(conj(flips))];		% Hyperecho refocusing



% -- These lines enable saving
global framenum filestem;
framenum=0;
filestem='/Users/brian/tmp/ims/frame';		% Place to save images!

Nspins = 36;

FZ = [0;0;1];	% Equilibrium

FZ = epg_animrf(FZ,pi/2,pi/2);		% 90 degree pulse.

for k=1:length(flips)
  FZ = epg_animgrad(FZ,Nspins);				% Crusher
  FZ = epg_animrf(FZ,abs(flips(k)),angle(flips(k)),Nspins);	% 90 deg pulse.
  FZ = epg_animgrad(FZ,Nspins);				% Crusher
  title('Spin Echo');
  pause(1);
  
    % -- If framenum and filestem exist, save to a file.
  for q=1:20	% Copy spin-echo frame many times to 'hold' in video.
    if (exist('epg_saveframe'))
      epg_saveframe;
    end;
  end;
end;




