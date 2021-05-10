% GET_STIM Summary of this function goes here
%
% Reference: "Peripheral Nerve Stimulation-Optimal GradientWaveform Design"
%            Rolf F. Schulte1 and Ralph Noeske
%            https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25440
%
%   Detailed explanation goes here...
%
% MLoecher@STANFORD.EDU (March 2020) for Rad229
% DBE@STANFORD.EDU (March 2021) for Rad229 - Minor updates

function stim_out = get_stim( G, dt )

%% Define the PNS model coefficients
L_grad = 0.333; % Effective gradient coil length [m] (sometimes called alpha)
r = 20.0; % Rheobase [T/s] {gradient coild specific}
% r = 23.4;
c = 334e-6; % Chronaxie time [s] {gradient coild specific}
Smin = r / L_grad; % Threshold at which 50% of people are stimulated

%%
coeff = zeros( numel(G) , 1);
for i = 1:(numel(G))
  coeff(i) = ( c / ((c + dt * ( numel(G) - 1 ) - dt * ( i - 1 ) ) ^ 2.0 ) / Smin );
end

stim_out = zeros(numel(G), 1);

for j = 1:(numel(G)-2)
  %         fprintf('%d \n', j);
  ss = 0;
  for i = 1:(j+1)
    %             fprintf('%d  %d  %d  %d  %d\n', j, (numel(G)-1), i, (j+1), numel(coeff)-j+i-1);
    %             fprintf('-----\n');
    ss = ss + coeff( numel(coeff)-j+i-1 ) * (G(i+1)-G(i));
  end
  stim_out(j) = ss;
end

end