function [FpFmZ] = epg_X_m0(f,opt)
%	Creates initial EPG-X magnetizations state.
%	
%	INPUT:
%		f = fraction of the second compartment
%       opt= put 'MT' to create the reduced MT case (4 components) or
%       empty for the full state
%
%	OUTPUT:
%		FpFmZ state.
%
    if nargin<2
        opt = 'FULL';
    end
    FpFmZ = zeros(6,1);
    FpFmZ(3)=(1-f);
    FpFmZ(6)=f;
    if strcmp(opt,'MT')
        FpFmZ(4:5)=[];
    end
end