%	function f = time2freq(t)
%
%	Function converts a time array to the corresponding frequency
%	points such that the (shifted) FFT of an array at times t leads
%	to frequency points f.  
%
%	Example:  time2freq([-5:4]) = [-.5:0.1:0.4]
%

function f = time2freq(t)

dt = t(2)-t(1);
df = 1/ (length(t))/dt;
np = length(t);

f = [-ceil((np-1)/2):floor((np-1)/2)]*df;

