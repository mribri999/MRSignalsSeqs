function [FpFmZ] = epg_X_grad(FpFmZ,noadd)
%	Propagate EPG-X states through a "unit" gradient.
%	
%	INPUT:
%		FpFmZ = 4xN or 6xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%
%	OUTPUT:
%		Updated FpFmZ state.
%
if (nargin < 2) noadd=0; end;	% Add by default.  
if size(FpFmZ,1) == 4 %MT case
    temp1 = epg_grad(FpFmZ(1:3,:),noadd);
    FpFmZ = [temp1; [FpFmZ(4,:) zeros(1,size(temp1,2)-size(FpFmZ,2))]];
else % full case
    temp1 = epg_grad(FpFmZ(1:3,:),noadd);
    temp2 = epg_grad(FpFmZ(4:6,:),noadd);
    FpFmZ = [temp1;temp2];
end

end

