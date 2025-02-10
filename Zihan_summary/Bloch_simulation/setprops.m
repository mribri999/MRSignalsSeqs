%
%	function setprops(handle, propertylist,debug)
%
%	Function sets the properties of the handle and children
%	according to propertylist.
%
%	propertylist is a list object of the form
%		{ type1, prop1, val1, type2, prop2, val2, ... }
%
%	For each child of handle, if the type matches typeN, and
%	there is a property matching propN, then it is set to valN.
%
%	B.Hargreaves
%

function setprops(handle, propertylist,debug)

% Do nothing... not working on 2016/2017 so commenting out!
