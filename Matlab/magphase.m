%
%	function magphase(x,arr)
%
%	Function plots in two subplots the magnitude and
%	phase of a complex-valued signal.
%

function magphase(x,arr)

if (nargin < 2)
	arr = x;
	x = 1:length(arr);
end;


mg = abs(arr);
ph = angle(arr);

subplot(2,1,1);
plot(x,mg);
ylabel('Magnitude');

subplot(2,1,2);
plot(x,ph/pi);
ylabel('Phase/\pi');
a = axis;
axis([a(1) a(2) -1 1]);
grid on;




