clear all
close all
%
Sim_times   = 1e2; % Monte-Carlo simulation iterations
%
SNR_dB      = -10:1:20; % in dB
SNR         = 10.^(SNR_dB./10);
%
RtF         = 1; % for User F
RtN         = 1; % for User N
R_on_off    = 1; % for ON/OFF
rho2        = 2^RtF - 1; % full time slot
rho1        = 2^RtN - 1;
rho2h       = 2^(2*RtF) - 1; % half time slot
rho1h       = 2^(2*RtN) - 1;
rho0        = 2.^R_on_off - 1; % SNR threshold for ON/OFF
% NOMA parameters
pF          = .8;
pN          = 1 - pF;
theta       = pF/pN;
% Path-loss model
L           = 1e3; % 30 dB path-loss at reference distance
epsilon     = 2.7; % path-loss exponent
d0          = 1; % reference distance
dSF         = 10;
dSN         = 3;
dNF         = dSF - dSN;
%
lSN         = L*((dSN/d0)^(-epsilon)); % lambda_SN
lNF         = L*((dNF/d0)^(-epsilon)); % lambda_NF
lSF         = L*((dSF/d0)^(-epsilon)); % lambda_SF
lrsi        = .1; % lambda_rsi
% Rician parameters
K_dB        = 3; % Rician K-factor in dB
K           = 10.^(K_dB./10);
% Rician Distribution parameters
nu          = sqrt(K/(K+1)*lrsi); % Noncentrality parameter
sigma       = sqrt(lrsi/2/(K+1)); % Scale parameter
%% SIMULATION
for ss = 1:length(SNR_dB)
    fprintf('SNR = %d dB \n',SNR_dB(ss))
    % Channel modeling
    hSF     = random('Rayleigh',sqrt(lSF/2),[1,Sim_times]);
    hSN     = random('Rayleigh',sqrt(lSN/2),[1,Sim_times]);
    hNF     = random('Rayleigh',sqrt(lNF/2),[1,Sim_times]);
    hrsi    = random('rician',nu,sigma,[1,Sim_times]); % Rician fading
    % Channel gains
    gSF     = abs(hSF).^2;
    gSN     = abs(hSN).^2;
    gNF     = abs(hNF).^2;
    grsi    = abs(hrsi).^2;
    % SNR/SINR Conventional
    snr_SF_co       = SNR(ss)*pF.*gSF ./ (SNR(ss)*pN.*gSF + 1);
    snr_SN_xF_co    = SNR(ss)*pF.*gSN ./ (SNR(ss)*pN.*gSN + 1);
    snr_SN_xN_co    = SNR(ss)*pN.*gSN;
    %
    OP_F_co_sim_temp = 0;
    OP_N_co_sim_temp = 0;
    %
    for zz = 1:Sim_times
        % for F
        if (snr_SF_co(zz) < rho2)
            OP_F_co_sim_temp = OP_F_co_sim_temp + 1;
        end
        % for N
        if (snr_SN_xF_co(zz) < rho2)
            OP_N_co_sim_temp = OP_N_co_sim_temp + 1;
        elseif (snr_SN_xF_co(zz) >= rho2) && (snr_SN_xN_co(zz) < rho1)
            OP_N_co_sim_temp = OP_N_co_sim_temp + 1;
        end
    end
    %
    OP_F_co_sim(ss) = OP_F_co_sim_temp/Sim_times;
    OP_N_co_sim(ss) = OP_N_co_sim_temp/Sim_times;
    %% Analysis
    a1      = SNR(ss)*pF;
    a2      = SNR(ss)*pN;
    a3      = SNR(ss);
    mu1     = rho1/a2;
    mu2     = rho2/(a1-a2*rho2);
    % User N
    if rho2 < theta
        if mu1 > mu2
            OP_N_co_ana(ss) = Psi(lSN,mu1);
        else
            OP_N_co_ana(ss) = Psi(lSN,mu2);
        end
    else
        OP_N_co_ana(ss) = 1;
    end
    % User F
    if rho2 < theta
        OP_F_co_ana(ss) = Psi(lSF,mu2);
    else
        OP_F_co_ana(ss) = 1;
    end
    %
end
%% Plot
semilogy(SNR_dB,OP_N_co_sim,'r>:',SNR_dB,OP_N_co_ana,'r-');
hold on
semilogy(SNR_dB,OP_F_co_sim,'bo:',SNR_dB,OP_F_co_ana,'b-');
legend('User N (sim.)', 'User N (ana.)', ...
    'User F (sim.)', 'User F (ana.)')
%
xlabel('SNR (dB)')
ylabel('Outage Probability (OP)')
%
save data_OP_N_co_sim.dat OP_N_co_sim -ascii
save data_OP_N_co_ana.dat OP_N_co_ana -ascii
save data_OP_F_co_sim.dat OP_F_co_sim -ascii
save data_OP_F_co_ana.dat OP_F_co_ana -ascii

