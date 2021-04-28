
%% Rad229_Motion_Artifacts_Demo – Demonstrates the impact of motion on both images and k-space.
%
%  This function is designed to demonstrate the effects of bulk motion and/or
%  pulsatile motion on MR image acquisition. The main feature is a digital
%  phantom with bulk and/or pulsatile motion options, plus noise. The
%  simulation is not a a precise MRI simulation per se owing to a simple 
%  treatment of noise and a lack of relaxation and general acqusition paramters.
%
% DBE@STANFORD.EDU (April 2020) for Rad229
% DBE@STANFORD.EDU (April 2021) for Rad229

%% Define MRI system constants
sys = Rad229_MRI_sys_config;

%% Set-up some default acquisition parameters
acq.Nx = 64;                            % Matrix is NxN (then padded later to accomodate motion)
acq.dt = 100e-3;                        % Time step [s]
acq.N_phs = 16;                         % Number of image phases [#]
acq.dur = acq.dt * ( acq.N_phs - 1 );   % Duration of motion [s]
acq.t = 0 : acq.dt : acq.dur;           % Time vector during motion [s]
acq.t_norm = ( 0 : acq.N_phs - 1 ) / acq.N_phs; % Normalized time vector

%% Generate the phantom object
[ P , M ] = Rad229_MRI_Phantom( acq );
P = repmat( P , [ 1 1 acq.N_phs ]);  % Replicate the object in time

%% Define a pulsatile "vessel" signal frequency and amplitude
obj.f_p = 2.0;   % Pulse frequency [cycles / phase] (e.g. 1 is one cycle per "beat")
obj.A_p = 1.1;   % Pulse amplitude [A.U.]
obj.num = 4;     % Identify the phantom sub-object (default phantom has 11 objects)

%% Create pulsating "vessel" object
obj.pulse_amp = obj.A_p * ( 1 + sin( 2 * pi * obj.f_p * acq.t_norm ) );
VSL = repmat( M( : , : , obj.num ) , [1 1 acq.N_phs] ); % Replicate object in time
VSL = VSL .* repmat( reshape( obj.pulse_amp , [ 1 1 length(obj.pulse_amp) ] ) , [acq.Nx acq.Nx 1] );

%% Define the bulk motion field for the object
obj.dx_max = round( 0 );     % Max x-displacement [pixels]
  warning('X-displacements not implemented.');
obj.dy_max = round( 8 );     % Max Y-displacement [pixels]
obj.f_bx = 0.5;   % Bulk frequency [cycles / phase] (e.g. 1 is one cycle per "beat")
obj.f_by = 0.5;   % Bulk frequency [cycles / phase] (e.g. 1 is one cycle per "beat")

% Define bulk motion path for each phase
obj.xt = obj.dx_max * ( sin( 2 * pi * obj.f_bx * acq.t_norm ) );
obj.yt = obj.dy_max * ( sin( 2 * pi * obj.f_by * acq.t_norm ) );

% Apply bulk motion and combine the pulsating "vessel" object with the static object
for n = 1 : acq.N_phs
  PHT( : , : , n) = circshift( P( : , : , n) + VSL( : , : , n) , round( obj.yt(n) ) , 1 );
end

%% Add noise (Note: this is an over-simplication of noise in MRI acquisition)
PHT = PHT + 0.05 * max(PHT(:)) * randn (size (PHT) );

%% Display the phantom motion
figure; imagescn(PHT);

PHT_scl = 255 * mat2gray( abs( PHT ) ) + 1; % Scale range for IMMOVIE function
P_mov = reshape( PHT_scl , [ size(PHT,1) size(PHT,2) 1 size(PHT,3) ] );
mov = immovie( P_mov , hot );
  implay(mov,4);

%% Compute the k-space for each phase
kspc = complex(zeros(size(PHT)));
for n=1:acq.N_phs
  kspc(:,:,n)=fftshift(fft2(PHT(:,:,n)));
end

figure; imagescn( abs( kspc ) ); colormap( hot );

%% Extract a the ky-lines acquired for each motion state (i.e. phase)
%  Because the object is moving we only acquire a few lines before it discretely moves again.
%  Note: The indexing indicates a top-down (row-by-row) acquistion.
VPS = size( P , 1 ) / acq.N_phs;  % Views-per-segment (lines of k-space acquired per phase)
for m = 1 : acq.N_phs
  ind0 = 1 + ( m - 1 ) * VPS;  % Index to first ky-line
  ind1 = m * VPS;              % Index to last ky-line
  kspc_motion( ind0 : ind1 , : ) = kspc( ind0 : ind1 , : , m ); % Copy lines from the mth phase
end

% Recover the acquired image for the moving/pulsatile object
PHT_motion=ifft2(fftshift(kspc_motion)); % Image of moving object

%% Generate figures...
set(groot,'defaultLineLineWidth',3,'defaultLineMarkerSize',25,'defaultAxesFontSize',15);
figure; hold on; colormap(hot);
subplot(2,2,1); hold on;
  plot(acq.t,obj.xt); plot(acq.t,obj.yt,'--'); plot(acq.t,obj.pulse_amp);
  title('Motion Path'); legend('X-path','Y-Path','Pulse','Location','NorthWest');

subplot(2,2,2); hold on;
  imagesc(abs(kspc_motion)); caxis([0 max(abs(kspc_motion(:)))/25]);
  title('Motion k-space MAG'); axis equal tight;

subplot(2,2,3); hold on;
  imagesc(abs(P(:,:,1))); caxis([0 1]);
  title('Object Mag.'); axis equal tight;
  
subplot(2,2,4); hold on;
  imagesc(abs(PHT_motion)); caxis([0 1]);
  title('Motion Object'); axis equal tight;
  
