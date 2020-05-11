%% This script demonstrates the impact of echo train length on sequence efficiency.
%
% Lecture_10_EPI_Sequence_Efficiency.m
%
% NOTE: Timing calculations are not precise for refocusing, pre-phasing,
% spoiling. They are meant to be approximate and representative.
%
% DBE@STANFORD.EDU (May 2020) for Rad229

%% Define some constants
gamma_bar=42.57e6;       % 1H gyromagnetic ratio [Hz/T]
dt=10e-6;                % Time-step        [s]
G_max = 100e-3;          % Gradient maximum [T/m]
S_max = 1500;             % Slewrate maximum [T/m/s]
dG=S_max*dt;             % Maximum gradient increment per time step

T2_star = 30e-3;         % T2-star for a typical intermediate tissue 
T1 = 500e-3;             % T1 for a typical intermediate tissue
FA = 20;                 % Flip angle for imaging

%% Define some acquisition parameters
TBW=32;                   % Time-bandwidth product
dur_RF=2e-3;             % RF pulse duration [s]
z_thick=10e-3;            % Slice thickness [m]

FOVx=0.256;              % Field-of-view [m]
Nx=128;                  % Number of x-pixels [#]

%% Design a slice-select gradient
Gss_max = (TBW/dur_RF)/(gamma_bar*z_thick);  % Gss plateau
Gss_rise = 0:dG:Gss_max;                     % Rise to Gss
Gss_fall = Gss_max:-dG:0;                    % Fall from Gss
Gss = [Gss_rise Gss_max*ones(1,numel(0:dt:dur_RF)-1) Gss_fall];     % Slice-select gradient
% Grf is a hack...should compare areas, but just used half the plateau duration
Grf = -[Gss_rise Gss_max*ones(1,numel(0:dt:dur_RF/2)-1) Gss_fall];  % Slice-select refocusing
Gss = [Gss Grf 0];                             % Composite gradient waveform

%% Design a readout gradient
rBW=16000;               % Receiver bandwidth [+/-Hz]
dt=1/(2*rBW);            % Sampling dwell time [s]
ETL=128;                 % Maximum echo train length (ETL)

Gro_max=(2*rBW)/(gamma_bar*FOVx);        % Readout gradient strength [T/m]
t_RO=Nx*dt;                              % Readout plateau duration [s]
t_PR=(Nx/2)*dt;                          % Pre-phaser plateau duration [s]
Gro_rise=0:dG:Gro_max;
Gro_fall=Gro_max:-dG:0;

Gpre=-[Gro_rise Gro_max*ones(1,numel(0:dt:t_PR)-1) Gro_fall];   % Readout pre-phaser
Gplat = [Gro_rise Gro_max*ones(1,numel(0:dt:t_RO)-1) Gro_fall]; % Readout gradient

%% Design a spoiler (somewhat arbitrary at 4x RO gradient)
Gsp_max=4*Gro_max;           % Spoiler gradient strength [T/m]
Gsp_rise=0:dG:Gsp_max;       % Spoiler gradient rise
Gsp_fall=[Gsp_max:-dG:0 0];  % Spoiler gradient fall
Gsp = [Gsp_rise Gsp_max*ones(1,numel(0:dt:t_RO/2)-1) Gsp_fall]; % Spoiler gradient

%% Assemble the readout gradient and capture some timing info
for e=1:ETL
  G_tmp=[];
  for f=1:e
    sgn = (-1)^(mod(f,2)+1);      % Plus/minus to toggle readout gradient
    G_tmp = [G_tmp sgn*Gplat];    % Append each readout gradient
    TE{e}(f) = dt*(numel(Gsp)+numel(G_tmp)-numel(Gplat)/2); % Calculate TE
  end
  Gro{e} = [0 Gpre G_tmp 0];                 % Final gradient waveform [T/m]
  TR(e) = dt*numel([Gss Gpre Gro{e} Gsp]);   % Sequence TR
  t_DAQ(e) = dt*numel(abs(Gro{e})==Gro_max); % Time spent acquiring data per TR
  eta(e) = t_DAQ(e)./TR(e);                  % Sequence efficiency
end

%% Assemble each sequence...
for e=1:ETL
  G_final{e}=[Gss Gro{e} Gsp];               % Composite gradient vector
  t_final{e}=0:dt:dt*(numel(G_final{e})-1);  % Time vector
end
    
%% Calculate each echo amplitude...
for e=1:ETL
  E1(e)=exp(-(TR(e)/T1));        % T1-relaxation
  E2{e}=exp(-(TE{e}/T2_star));   % T2-star relaxation
  % Steady-state gradient echo signal equation
  A_echo{e} = ((1-E1(e))/(1-cosd(FA)*E1(e)))*(sind(FA)*E2{e});
end

%% Figure of some multi-echo composite gradients
set(groot,'defaultLineLineWidth',3,'defaultLineMarkerSize',25,'defaultAxesFontSize',15);
figure; hold on;
  subplot(3,1,1); hold on; axis([0 0.08 -0.05 0.05]);
    plot(t_final{4},G_final{4});  
    plot(0:dt:dt*(numel(Gss)-1),Gss);
    ylabel('Gradient [T/m]'); title('Multi-Echo Sequence');
    legend('Readout Gradient','Slice Selection');
  subplot(3,1,2); hold on; axis([0 0.08 -0.05 0.05]);
    plot(t_final{8},G_final{8});   
    plot(0:dt:dt*(numel(Gss)-1),Gss);
    ylabel('Gradient [T/m]');
  subplot(3,1,3); hold on; axis([0 0.08 -0.05 0.05]);
    plot(t_final{16},G_final{16}); 
    plot(0:dt:dt*(numel(Gss)-1),Gss);
    xlabel('Time [s]'); ylabel('Gradient [T/m]');

%% Figure of Multi-Echo Sequence Performance
ETL0=32; ETL1=ETL;
f=figure; hold on; title('Multi-Echo Sequence Performance'); 
yyaxis left; plot(1:ETL,eta);
  ylabel('Sequence Efficiency [\eta]');
yyaxis right; plot(1:ETL,A_echo{ETL}); plot(1:16,A_echo{16});
  ylabel('Signal Amplitude [A.U.]');
  xlabel('Echo Train Length [#]'); 
  legend('Efficiency [\eta]',['A_{Echo} (ETL=',num2str(ETL1),')'],['A_{Echo} (ETL=',num2str(ETL0),')'],'Location','east');