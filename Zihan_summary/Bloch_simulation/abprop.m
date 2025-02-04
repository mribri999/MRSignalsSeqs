%
%	function [A,B,mss] = abprop(A1,B1,A2,B2,A3,B3,...)
%
%	Function 'propagates' matrices A and B to give overall A,B
%	such that m' = A*m+B.  The Ai,Bi are applied in order, ie A1,B1
%	are the first applied to m, then A2,B2 and so on.
%
%	If mss is provided, the steady-state is calculated - in which case
%	the A1,B1 list should start right after the point where it is desired
%	to calculate mss.
%
%	If an Ai is 3x4, then it is assumed to be [Ai Bi]
%
%	If a Bi vector is omitted (the next argument is 3x3 or 3x4, 
%	it is assumed to be zero.
%
function [A,B,mss] = abprop(varargin)

A=eye(3); B=[0;0;0];		% Start with identity
mss = [];			% Default

argcount=1;
while (argcount <= nargin)
  Ai = varargin{argcount};
  argcount = argcount+1;

  if (size(Ai) == [3 3])	% If next argument is 3x1, that's Bi,
				% otherwise (not 3x1 or end) Bi=0.
    if (argcount <= nargin) 
      Bi = varargin{argcount};
      if (size(Bi)==[3 1]) 
	argcount = argcount+1;
      else
	Bi = [0;0;0];	% Wrong size for Bi, so ignore next argument
      end;
    else
      Bi = [0;0;0];
    end;
  end;

  if (size(Ai) == [3 4])	% -- If 3x4, then of form [Ai Bi]
    Bi = Ai(:,4);
    Ai = Ai(:,1:3);
  end;

  % -- If Ai is not 3x4 or 3x3 then force an error.
  if (size(Ai,1)~=3) || (size(Ai,2)~=3)	
    tt = sprintf('Argument %d should be 3x3 or 3x4',2*count-1);
    error(tt); 
  end;

  %disp('Propagating...');	% -- Display for testing
  %disp ([ Ai Bi]); 

  A = Ai*A;
  B = Ai*B+Bi;

end;  

if (nargout>2) 			% -- Try to calculate Steady-State Mag.
  if (rank(eye(3)-A)<3)
    warning('No Unique Steady-State - mss will be []');
  else
    mss = inv(eye(3)-A)*B;
  end;
end;


