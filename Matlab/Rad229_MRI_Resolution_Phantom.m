%% Rad229_MRI_Resolution_Phantom – Create a resolution phantom with some noise.
%
% SYNTAX  - [ P ] = Rad229_MRI_Resolution_Phantom( acq , obj )
%
% INPUTS  - acq is a structure that needs to contain acq.Nx (number of
%             pixels defining the *square* phantom matrix size.
%
%         - obj is a structure that needs to contain a obj.obj.T2star vector
%             of FIVE T2-star values in seconds.
%
% OUTPUTS - P - Is the phantom object matrix [acq.Nx x acq.Nx] with assigned T2-star values
%
% DBE@STANFORD.EDU (April 2021) for Rad229

function P = Rad229_MRI_Resolution_Phantom( acq , obj )

if nargin == 0
  acq.Nx = 128; % Matrix will be NxN
  obj.T2star = [10e-3 25e-3 50e-3 100e-3 1000e-3]; % Range of T2-star values [s]
end

P = zeros(acq.Nx);
R0 = 16 :  2 : acq.Nx  
R1 = 16 :  2 : acq.Nx % Fine line spacing
C0 =  9 : 24 : acq.Nx  
C1 = 25 : 24 : acq.Nx % 
for c = 1 : 5
  for r = 1 : 5
    P( R0(r) : R1(r) , C0(c) : C1(c) ) = obj.T2star(r); 
  end 
end

R0 = 48 :  4 : acq.Nx;  R1 = 49 :  4 : acq.Nx; % Medium line spacing
C0 =  9 : 24 : acq.Nx;  C1 = 25 : 24 : acq.Nx;
for c = 1 : 5
  for r = 1 : 5
    P( R0(r) : R1(r) , C0(c) : C1(c)) = obj.T2star(r); 
  end
end

R0 = 88 :  8 : acq.Nx;  R1 = 91 :  8 : acq.Nx; % Coarse line spacing
C0 =  9 : 24 : acq.Nx;  C1 = 25 : 24 : acq.Nx;
for c=1:5
  for r=1:5
    P( R0(r) : R1(r) , C0(c) : C1(c) ) = obj.T2star(r); 
  end
end

% P = zeros(acq.Nx);
% C0 = 16 :  2 : acq.Nx;  C1 = 16 :  2 : acq.Nx; % Fine line spacing
% R0 =  9 : 24 : acq.Nx;  R1 = 25 : 24 : acq.Nx; % 
% for c = 1 : 5
%   for r = 1 : 5
%     P( R0(r) : R1(r) , C0(c) : C1(c) ) = obj.T2star(r); 
%   end 
% end
% 
% C0 = 48 :  4 : acq.Nx;  C1 = 49 :  4 : acq.Nx; % Medium line spacing
% R0 =  9 : 24 : acq.Nx;  R1 = 25 : 24 : acq.Nx;
% for c = 1 : 5
%   for r = 1 : 5
%     P( R0(r) : R1(r) , C0(c) : C1(c)) = obj.T2star(r); 
%   end
% end
% 
% C0 = 88 :  8 : acq.Nx;  C1 = 91 :  8 : acq.Nx; % Coarse line spacing
% R0 =  9 : 24 : acq.Nx;  R1 = 25 : 24 : acq.Nx;
% for c=1:5
%   for r=1:5
%     P( R0(r) : R1(r) , C0(c) : C1(c) ) = obj.T2star(r); 
%   end
% end
% 
% P=P'; % Want the resolution elements parallel to phase encode direction

return