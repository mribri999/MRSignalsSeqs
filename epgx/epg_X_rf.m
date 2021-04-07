function [FpFmZ,RR] = epg_X_rf(FpFmZ,alpha,phi,WT)
%	Simulates action of RF pulse on EPG-X states
%	
%	INPUT:
%		FpFmZ = 4xN or 6xN vector of F+, F- and Z states.
%       alpha: flip angle in radians
%       phi: phase of the RF pulse (radians)
%       WT: saturation of the bound pool componen will be exp(-WT), ignore for
%       full case
%	OUTPUT:
%		FpFmZ state.
%       RR: rotation matrix

    if size(FpFmZ,1) == 4 %MT case      
        [temp1,RR] = epg_rf(FpFmZ(1:3,:),alpha,phi);
        FpFmZ = [temp1; exp(-WT)*FpFmZ(4,:)];
    else % full case
        [temp1,RR] = epg_rf(FpFmZ(1:3,:),alpha,phi);
        [temp2,RR] = epg_rf(FpFmZ(4:6,:),alpha,phi);
        FpFmZ = [temp1;temp2];
    end
end

