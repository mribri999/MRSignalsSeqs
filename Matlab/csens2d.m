%
%	function [ncs,cs] = csens2d(klow)
%
%	Function calculates 2D coil sensitivities from an array klow,
%	that is a calibration region (ie mostly zero except near the 
%	center.
%
%	INPUT:
%		klow = Nx x Ny x Nc matrix
%
%	OUTPUT:
%		ncs = coil sensitivities, normalized by RMS
%		cs = coil sensitivities, raw
%		cmrms = RMS of sensitivities


function [ncs,cs,cmrms] = csens2d(klow)

sz = size(klow);
nc = sz(3);

cs = fftshift(fft(fftshift(klow,1),[],1),1);    % FFT x
cs = fftshift(fft(fftshift(cs,2),[],2),2);      % FFT y
cmrms = sqrt( sum( conj(cs).*cs , 3) );         % RMS of maps, coil dimension
meancm = mean(cmrms(:));

% -- Normalize
% -- If RMS is 0, replace values in channel 1 
f = find(cmrms(:)==0);
cs = reshape(cs,sz(1)*sz(2),nc);		% Npix x Nc
cmrms(f)=meancm/100;

cs(f,:) = meancm/10000;			% Small values in background
ncs = cs ./ (cmrms(:)*ones(1,nc));	% Normalized matrix
ncs = reshape(ncs,sz);
cs = reshape(cs,sz);


 
