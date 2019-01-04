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
R_on_off    = 2; % for ON/OFF
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
    pout_F_s3_sim_temp = 0;
    %
    count1 = (snr_SN_xF_co < rho2h) & (snr_SF_co < rho2h);
    count2 = (snr_SN_xF_co >= rho2h) & ...
        (snr_SF_co + snr_NF_HD < rho2h);
    OP_F_s3_sim(ss) = sum(count1+count2)/Sim_times;
    %% ANALYSIS
    a1 = SNR(ss)*pF;
    a2 = SNR(ss)*pN;
    a3 = SNR(ss);
    chi = 1/(a2*a3*lNF*lSF);
    xi = a1/a2;
    mu2h = rho2h/(a1-a2*rho2h);
    %
    Upsilon_mu2h_rho2h = 1 - exp(-rho2h/lSF/(a1-a2*rho2h))  ...
        - 1/a2/a3/lNF/lSF*exp(-rho2h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*rho2h/(a1-a2*rho2h)+a3*lNF,chi,xi));
    %
    if rho2h < theta
        OP_F_s3_ana(ss) = Psi(lSN,mu2h)*Psi(lSF,mu2h) ...
            + (1-Psi(lSN,mu2h))*Upsilon_mu2h_rho2h;
    else
        OP_F_s3_ana(ss) = 1;
    end
end
%
semilogy(SNR_dB,OP_F_s3_sim,'s');
hold on
semilogy(SNR_dB,OP_F_s3_ana,'*-');
% axis([-10 30 1e-6 1e0])
grid on
%
legend('User F, Scheme III (sim.)', 'User F, Scheme III (ana.)')
xlabel('SNR (dB)')
ylabel('Outage Probability')
%
save data_OP_F_s3_sim.dat OP_F_s3_sim -ascii
save data_OP_F_s3_ana.dat OP_F_s3_ana -ascii