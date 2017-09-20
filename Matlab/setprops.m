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

if (nargin < 1)
	handle = gcf;
end;
if (nargin < 2)
	propertylist = bhprops;
end;
if (nargin < 3)
	debug = 0;
end;
proptot = 0;

propertylist = propertylist(:);

h = printchildren(handle,0,1,1);	% Last 1->0 to print.

lp = floor(length(propertylist)/3);

for k=1:length(h)
	hh = h(k);
	htype = get(hh,'Type');
	for q=1:lp
		ptype = propertylist{3*q-2};
		pprop = propertylist{3*q-1};
		pval  = propertylist{3*q  };
		if (debug >1)
			{ptype pprop pval}
		end;
		if (length(htype)==length(ptype)) % lengths match
		    match=1;
		    for p=1:length(htype)
			if (htype(p)~=ptype(p))
				match=0;
			end;
		    end;
		    if (match==1)
		        try
			    found=1;
			    if (debug)
			        tt=sprintf('Trying to Set %s for object %f (%s)',pprop, hh, htype);
				disp(tt);
			    end;
			    set(hh,pprop, pval);
		        catch
			    found=0;
		        end;	
			if (found)
			    proptot = proptot + 1;
		            if (debug)
			        tt=sprintf('Set %s for object %f',pprop, hh);
			        disp(tt);
		    	    end;
		        end;
		    end;
		end;
	end;
end;
if (debug)
  tt = sprintf('Found %d objects. ',length(h));
  disp(tt);
  tt = sprintf('Set %d properties.',proptot);
  disp(tt);
end;


% --------------------------------------

function props = bhprops()

%	Returns a list of properties for setprops.

props = {
%	Type	Property    Value
	'axes', 'FontName', 'Helvetica',	% Axis numbering font
	'axes', 'FontSize', 12,		% Axis numbering size
	'text', 'FontName', 'Helvetica',	% Label/Title font
	'text', 'FontSize', 16',	% Label/Title font
	'line', 'LineWidth', 2 ,
	'axes', 'LineWidth', 1 
	};

props = props.';
props = props(:);




%	function hc=printchildren(handle, indent,num, noprint)
%
%	Function starts with a handle and tries to pull out all
%	the children of that handle, and their children, and so on.
%
%	Prints a list of the objects, types.  
%
%	Returns a single vector of the handles.


function hc=printchildren(handle, indent,num, noprint)
%	Prints the "type" of each child, and their children,
%	recursively.

if (nargin < 2)
	indent = 0;
end;
if (nargin < 3)
	num = 1;
end;
if (nargin < 4)
	noprint=0;
end;

deltaindent = 4;

%  Make a string of spaces with correct indent.
indentspace = sprintf('%3d) ',num);
for k=1:indent
	indentspace = sprintf('%s ',indentspace);
end;

%  Get type and children for handle, and print.
htype = get(handle,'Type');
tt = sprintf('%s%f -- %s',indentspace,handle,htype);
if (noprint==0)
	disp(tt);
end;
h1 = get(handle,'Children');

% The following because Matlab is evil, making some things
%	act like children, but be linked via other fields...
try h2=get(handle,'Title');
	catch h2=[];
end;
try h3=get(handle,'XLabel');
	catch h3=[];
end;
try h4=get(handle,'YLabel');
	catch h4=[];
end;
try h5=get(handle,'ZLabel');
	catch h5=[];
end;
h = [h1;h2;h3;h4;h5];

%  Recursively Print types for children.
hc = handle;

if (length(h)>0)
	for q=1:length(h);
		hc = [hc; printchildren(h(q),indent+deltaindent,num+length(hc),noprint)];
	end;
	if (noprint==0)
		disp(' ');	% Blank line after last child.
	end;
end;


