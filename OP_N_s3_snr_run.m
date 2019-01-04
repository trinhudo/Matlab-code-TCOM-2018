clear all
close all
%
Sim_times   = 1e3; % Monte-Carlo simulation iterations
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
epsilon     = 2.7;  % path-loss exponent
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
chi         = sqrt(K/(K+1)*lrsi); % Noncentrality parameter
sigma       = sqrt(lrsi/2/(K+1)); % Scale parameter
%% SIMULATION
for ss = 1:length(SNR_dB)
    fprintf('SNR = %d dB \n',SNR_dB(ss))
    % Channel modeling
    hSF     = random('Rayleigh',sqrt(lSF/2),[1,Sim_times]);
    hSN     = random('Rayleigh',sqrt(lSN/2),[1,Sim_times]);
    hNF     = random('Rayleigh',sqrt(lNF/2),[1,Sim_times]);
    % Channel gains
    gSF     = abs(hSF).^2;
    gSN     = abs(hSN).^2;
    gNF     = abs(hNF).^2;
    %
    % SNR/SINR Conventional
    snr_SF_co       = SNR(ss)*pF.*gSF ./ (SNR(ss)*pN.*gSF + 1);
    snr_SN_xF_co    = SNR(ss)*pF.*gSN ./ (SNR(ss)*pN.*gSN + 1);
    snr_SN_xN_co    = SNR(ss)*pN.*gSN;
    %
    snr_NF_HD      = SNR(ss).*gNF ;
    snr_SNF_HD     = min(snr_SN_xF_co,snr_SF_co+snr_NF_HD); % using MRC
    %
    OP_N_s3_sim_temp = 0;
    % Transition probabilites
    for zz = 1:Sim_times
        % for N
        if (snr_SN_xF_co(zz) < rho2h)
            OP_N_s3_sim_temp = OP_N_s3_sim_temp + 1;
        elseif (snr_SN_xF_co(zz) >= rho2h) && (snr_SN_xN_co(zz) < rho1h)
            OP_N_s3_sim_temp = OP_N_s3_sim_temp + 1;
        end
    end
    OP_N_s3_sim(ss) = OP_N_s3_sim_temp/Sim_times;
    %% ANALYSIS
    a1 = SNR(ss)*pF;
    a2 = SNR(ss)*pN;
    a3 = SNR(ss);
    %
    if rho2h < theta
        N31 = Psi(lSN,rho2h/(a1-a2*rho2h));
    else
        N31 = 1;
    end
    %
    if (rho2h < theta) && (rho2h/(a1-a2*rho2h) < rho1h/a2)
        N32 = Psi(lSN,rho1h/a2) - Psi(lSN,rho2h/(a1-a2*rho2h));
    else
        N32 = 0;
    end
    %
    OP_N_s3_ana(ss) = N31 + N32;
end
%
semilogy(SNR_dB,OP_N_s3_sim,'o:');
hold on
semilogy(SNR_dB,OP_N_s3_ana,'-');
%
legend('User N, Scheme III (sim.)')
title('OP of User N in Scheme III')
xlabel('SNR (dB)')
ylabel('Outage Probability')

save data_OP_N_s3_sim.dat OP_N_s3_sim -ascii
save data_OP_N_s3_ana.dat OP_N_s3_ana -ascii