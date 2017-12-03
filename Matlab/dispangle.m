%
%	function dispangle(arr)
%
%	Function displays the angle of the complex data
%	on a scale of black=-pi to white=pi.
%

function dispangle(arr)

angarr = angle(arr)+pi;
dispim(angarr,0,2*pi);

