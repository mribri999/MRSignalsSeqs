%% Rad229_MRI_Phantom – Create an object with MRI related properties.
%
% SYNTAX  - [ P , M ] = Rad229_MRI_Phantom(acq)
%
% INPUTS  - acq is a structure that needs to contain acq.Nx (number of
%           pixels defining the *square* phantom matrix size.
%
% OUTPUTS - P - Is the phantom object matrix [acq.Nx x acq.Nx]
%           M - A 3D logical mask matrix. Each layer is an object in the phantom.
%
% DBE@STANFORD.EDU (April 2021) for Rad229

function [ P , M ] = Rad229_MRI_Phantom(acq)

if nargin == 0
  acq.Nx = 128;                           % Matrix is NxN (then padded later to accomodate motion)
end

%%  Use the default "modified shepp-logan" phantom parameters and define "tissue" parameters
[ P , E ] = phantom( 'modified shepp-logan' , acq.Nx ); % Create a default phantom

%% Create a separate object for each feature
for n = 1 : size( E , 1 )
  tmp = phantom( [ 1 E( n , 2 : end ) ] , acq.Nx );  % Force each mask intensity to be ONE
  M( : , : , n ) = tmp ~= 0; % Logical mask for each object in phantom (objects may overlap)
end

% Parse the background and define some "white-matter" (bulk of phantom) and a "skull" (outer ring of "tissue")  
Q( : , : , 1 ) = ~M( : , : , 1 ) & ~M( : , : , 2 );  % Background
Q( : , : , 2 ) =  M( : , : , 1 ) &  M( : , : , 2 );  % White-matter
Q( : , : , 3 ) =  M( : , : , 1 ) & ~M( : , : , 2 );  % "Skull" or "fat"
M = cat( 3 , Q , M( : , : , 3 : end ) );
M( : , : , 2 ) =~ ( sum( M , 3 ) > 1 ) & Q( : , : , 2 );  % This makes the "white matter" only everywhere around the ROIs (no overlap)

return

%% NOT IMPLEMENTED YET...
% obj.Nx = 
% obj.pd(n).
% obj.T1(n).
% obj.T2(n).
% obj.fat(n).
% obj.off(n).
% obj.T2s(n).
% obj.obj(n) = 
% obj.tissue(n) = 

%% Assign tissue parameters to each feature
% T1_fat=260;
% T2_fat=85;
% PD_fat=1.0;
% 
% T1_WM=790;
% T2_WM=92;
% PD_WM=0.9;
% 
% T1_CSF=2400;
% T2_CSF=180;
% PD_CSF=1.0;
% 
% T1_HEM=400;
% T2_HEM=50;
% PD_HEM=0.75;