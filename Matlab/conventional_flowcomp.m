function [FC, M0S, M1S, t_ss, G_ss, M0_f, M1_f] = conventional_flowcomp(params)

% Gss_plat = linspace( params.g_ss * 1e-3 , params.g_ss * 1e-3, ceil( params.p_ss * 1e-3 / params.dt ) );
Gss_plat = linspace( params.g_ss , params.g_ss , ceil( params.p_ss / params.dt ) );
Gss_down = linspace( params.g_ss , params.g_ss / 10 , 10 );
warning('Ramp to zero?')

G_ss = [ Gss_plat Gss_down ];

t_ss = length(G_ss) * params.dt;
t_vec = 0 : params.dt : ( length(G_ss) - 1) * params.dt;

% Calculate the gradient moments
% M0S = cumsum( ( G_ss .* params.dt ) .* 1e6 , 2);           % [mT/m x ms]
% M1S = cumsum( ( G_ss .* t_vec * params.dt) .* 1e9 , 2);    % [mT/m x ms^2]
M0S = cumsum( ( G_ss .* params.dt ) , 2);           % [T/m x s]
M1S = cumsum( ( G_ss .* t_vec * params.dt) , 2);    % [T/m x s^2]

% Store the last gradient moment value
M0S_ = M0S( 1 , end );
M1S_ = M1S( 1 , end );

% Design the flow comp gradient waveform
% for r = 1e-03 : 1e-03 : 5
%   r_ = ( ceil( r / params.dt / 1000 ) );
% for r = 1 : 1 : 5e3
for r = 1e-06 : 1e-06 : 5e-03
r_ = ( ceil( r / params.dt ) );
  h = r * params.smax;
  M02 = ( -h * r + sqrt( ( h * r ) ^ 2 + 2 * ( h * r * M0S_ + M0S_^2 + 2 * h * M1S_ ) ) ) /2;
  M01 = M02 + M0S_;
  w1 = M01 / h + r;
  w2 = M02 / h + r;
%   w1_ = ceil( w1 / params.dt / 1000 ) ;
%   w2_ = ceil( w2 / params.dt / 1000 ) ;
  w1_ = ceil( w1 / params.dt ) ;
  w2_ = ceil( w2 / params.dt ) ;
  
  if  (w1_ - 2 * r_ <= 1) || ( w2_ - 2 * r_ <= 1)
    h1 = M01 / ( r_ - w1_ ) * 100;
    h2 = -M02 / ( r_ - w2_ ) * 100;
    break
  end
end

% Construct the gradient waveform segments
Gfc_down1 = linspace( 0 , h1 , r_ );
Gfc_plat1 = linspace( h1 , h1 , w1_ - 2 * r_ );
Gfc_up = linspace( h1 , h2 , 2 * r_ );
Gfc_plat2 = linspace( h2 , h2 , w2_ - 2 * r_ );
Gfc_down2 = linspace(h2,0,r_);

% Concatenate the gradient waveform segments
G = [Gfc_down1 Gfc_plat1 Gfc_up Gfc_plat2 Gfc_down2];

% FC = horzcat( G_ss * 1000 , G );

FC = horzcat( G_ss , G );

% Compute final gradient moments
t_vec = 0 : params.dt : ( length( FC ) - 1 ) * params.dt;
% M0_f = cumsum( ( FC .* params.dt ) .* 1e6 , 2 );           % [mT/m x ms]
% M1_f = cumsum( ( FC .* t_vec * params.dt ) .* 1e9 , 2);    % [mT/m x ms^2]
M0_f = cumsum( ( FC .* params.dt ) , 2 );           % [T/m x s]
M1_f = cumsum( ( FC .* t_vec * params.dt ) , 2);    % [T/m x s^2]

return