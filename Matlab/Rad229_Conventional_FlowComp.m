%% Rad229_Conventional_FlowComp – Demonstrate designing a flow-compensated gradient waveform.
%
%  Based, in part, on this paper:
%
%    Encoding strategies for three-direction phasecontrast MR imaging of flow. 
%    Pelc NJ, Bernstein MA, Shimakawa A, Glover GH. J Magn Reson Imaging 1991;1(4):405-413.
% 
%  SYNTAX - [G_FC, M0S, M1S, t_ss, G_ss, M0_FC, M1_FC] = Rad229_Conventional_FlowComp(params)
%
%  INPUTS -  params
%
%  OUTPUTS – G_FC
%            M0S
%            M1S
%            t_ss
%            G_ss
%            M0_FC
%            M1_FC
%            
% MLoecher@STANFORD.EDU (March 2020) for Rad229
% DBE@STANFORD.EDU (March 2021) for Rad229 - Minor updates

function [G_FC, M0S, M1S, t_ss, G_ss, M0_FC, M1_FC, t_FC] = Rad229_Conventional_FlowComp(params)

% Build the slice select gradient
Gss_plat = linspace( params.g_ss , params.g_ss , ceil( params.p_ss / params.dt ) );
% Gss_down = linspace( params.g_ss , params.g_ss / 10 , 10 );
Gss_down = linspace( params.g_ss , 0 , 10 );

G_ss = [ Gss_plat Gss_down ];

t_ss = length(G_ss) * params.dt;
t_vec = 0 : params.dt : (length(G_ss)-1)*params.dt;

% Get moments of the slice select that we need to compensate
M0S = cumsum( (G_ss .* params.dt) , 2);            % [T/m x s]
M1S = cumsum( (G_ss .* t_vec * params.dt) , 2);    % [T/m x s^2]

% Moments at end of slice select
M0S_ = M0S( 1 , end );
M1S_ = M1S( 1 , end );

% Now solve for the next two lobes as presented in:
%   Encoding strategies for three-direction phasecontrast MR imaging of flow. 
%   Pelc NJ, Bernstein MA, Shimakawa A, Glover GH. J Magn Reson Imaging 1991;1(4):405-413.
% But we add an additional search over risetime (r) to allow for waveforms with triangle lobes.

% Maximum riseitme [s]
rmax = params.gmax/params.smax;

% r = risetime in [s]
for r = 1e-06 : 1e-06 : rmax
  r_ = (ceil(r/params.dt)); % risetime in 'dt' units
  h = r*params.smax; % gradient strength [T/m]
  
  M02 = (-h*r+sqrt((h*r)^2+2*(h*r*M0S_ + M0S_^2 + 2*h*M1S_)))/2; % Area of the second lobe (see above paper for derivation) [T/m x s]
  M01 = M02 + M0S_; % Area of the first lobe (see above paper for derivation) [T/m x s]
  
  w1 = M01/h + r; % Duration of first lobe [s]
  w2 = M02/h + r; % Duration of second lobe [s]
  
  w1_ = (ceil(w1/params.dt)); % Duration of first lobe in 'dt' units
  w2_ = (ceil(w2/params.dt)); % Duration of first lobe in 'dt' units
  
  % Signed amplitudes [T/m]
  h1 = -h;
  h2 = h;
  
  % We stop this loop if the flattop time is less than one 'dt' unit (triangle lobe)
  if  (w1_-2*r_ <= 1) || (w2_-2*r_ <= 1)
    h1 = M01/(r_-w1_)*100;
    h2 = -M02/(r_-w2_)*100;
    break
  end
end

% Build the flow compensated bipolar
G = horzcat(linspace(0,h1,r_),linspace(h1,h1,w1_-2*r_),linspace(h1,h2,2*r_),linspace(h2,h2,w2_-2*r_),linspace(h2,0,r_));

% Output waveform in [T/m]
G_FC = horzcat(G_ss,G);

% Compute final gradient moments
t_vec = 0 : params.dt : ( length( G_FC ) - 1 ) * params.dt;
M0_FC = cumsum( ( G_FC .* params.dt ) , 2 );           % [T/m x s]
M1_FC = cumsum( ( G_FC .* t_vec * params.dt ) , 2);    % [T/m x s^2]

% Define a time vector
t_FC = ( 0 : (length(G_FC)-1) ) * params.dt;
end