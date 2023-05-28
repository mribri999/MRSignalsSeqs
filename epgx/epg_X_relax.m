function [FpFmZ] = epg_X_relax(FpFmZ,T,T1,T2,ka,deltab,f)
%	Relaxation and exchange of  EPG-X.
%	
%	INPUT:
%		FpFmZ = 4xN or 6xN vector of F+, F- and Z states.
%       T: time for the relaxation
%       T1: 1x2 vector with the T1 constant of both comparment
%       T2: 1x2 vector or scalar with the T2 constants
%       ka: forward exchange rate
%       deltab: frequency offset of 2nd comparment
%       f: fraction of second compartment
%
%	OUTPUT:
%		FpFmZ state.
%
   M0b = f; 
   M0a = 1-f; 
   kb = ka * M0a/M0b; % conservation of magnetization
 
   R1a=1/T1(1);
   R1b=1/T1(2); 
   R2a=1/T2(1);
   LambdaL = expm([-R1a-ka kb; ka -R1b-kb]*T); %longitudinal relaxation and exchange
   E1 = exp(-T./T1');
   Zoff = [M0a;M0b].*(1-E1); % longitudinal recovery 
   if size(FpFmZ,1) == 4 %MT case
        E2a = exp(-T*R2a);
        FpFmZ(1:2,:) = diag([E2a E2a])*FpFmZ(1:2,:);
        FpFmZ(3:4,:) = LambdaL*FpFmZ(3:4,:);
        FpFmZ(3:4,1) = FpFmZ(3:4,1)+Zoff;
   else % full case
        R2b=1/T2(2); 
        LambdaT = expm(T*[-R2a-ka 0 kb 0; 0 -R2a-ka 0 kb; ka 0 -R2b-kb-2*pi*1i*deltab 0; 0 ka 0 -R2b-kb+2*pi*1i*deltab]); %xversal relaxation and exchange
        temp1 = LambdaT*FpFmZ([1;2;4;5],:); 
        temp2 = LambdaL*FpFmZ([3;6],:);
        temp2(:,1) = temp2(:,1) + Zoff;
        FpFmZ = [temp1(1:2,:);temp2(1,:);temp1(3:4,:);temp2(2,:)];
   end

end

