%
%	Function b = bvalue(gradwave,T)
%
%	Function calculates the b value of a gradient waveform.
%	If there is a 180 pulse, then flip the gradient sign on one
%	side or another!
%
%	Input:
%		gradwave = waveform, in mT/m
%		T = sampling period in seconds
%
%	Output:
%		b = b value in s/mm^2
%

function b = bvalue(gradwave,T)

gamma = 2*pi*42.58;		% rad * Hz/T
intg = cumsum(gradwave)*T;	% mT/m*s
b = gamma^2 * sum( intg.^2) * T;	% rad^2 /T^2 * (mT^2 / m^2) s
				


