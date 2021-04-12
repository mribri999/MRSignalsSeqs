% Rad229_Fourier_Encoding_Demo.m is a script and directions to help demonstrates MRI Fourier 
% encoding (k-space).
%
% DBE@STANFORD.EDU (March 2021) for Rad229
error('This function is incomplete.')
%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Load some MRI data (default with MATLAB installation)
load MRI; D = double(D(:,:,1,10)); % Grab a single slice and convert to DOUBLE
I = imresize(D, 0.3);              % It will be helpful to use a low-res object
I = I(:,2:end-1);                  % Assymetric matrix dimensions are helpful for debugging!


%% Define the acquistion parameters that are well matched to the object "I"
%  acq.FOVx = ??? etc.


%% Use Rad229_Fourier_Encoding.m (or similar) to compute all of the Fourier sampling 
%  coefficients (F) for the object "I" using "acq"


%% Compute each k-point for each Fourier encoding pattern and store in a matrix


%% Compute the inverse FFT of your k-space data to recover an image of the object


