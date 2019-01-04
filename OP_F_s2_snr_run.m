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
L           = 1e3; % -30 dB path-loss at reference distance
epsilon     = 2.7; % path-loss exponent
d0          = 1; % reference distance
dSF         = 10;
dSN         = 3;
dNF         = dSF - dSN;
%
lSN         = L*((dSN/d0)^(-epsilon)); % lambda_SN
lNF         = L*((dNF/d0)^(-epsilon)); % lambda_NF
lSF         = L*((dSF/d0)^(-epsilon)); % lambda_SF
% Rician parameters
K_dB        = 25; % Rician K-factor in dB
K           = 10.^(K_dB./10);
lrsi        = K/(K+1); % lambda_rsi
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
    hrsi    = random('rician',chi,sigma,[1,Sim_times]); % Rician fading
    % Channel gains
    gSF     = abs(hSF).^2;
    gSN     = abs(hSN).^2;
    gNF     = abs(hNF).^2;
    grsi    = abs(hrsi).^2;
    % SNR/SINR
    snr_SF_co       = SNR(ss)*pF.*gSF ./ (SNR(ss)*pN.*gSF + 1);
    snr_SN_xF_co    = SNR(ss)*pF.*gSN ./ (SNR(ss)*pN.*gSN + 1);
    snr_SN_xN_co    = SNR(ss)*pN.*gSN;
    %
    snr_SN_xF_FD   = SNR(ss)*pF.*gSN ./ (SNR(ss)*pN.*gSN + SNR(ss).*grsi + 1);
    snr_SN_xN_FD   = SNR(ss)*pN.*gSN ./ (SNR(ss).*grsi + 1);
    snr_NF_FD      = SNR(ss).*gNF ./ (SNR(ss).*gSF + 1);
    snr_SNF_FD     = min(snr_SN_xF_FD,snr_NF_FD);
    % Transition probabilites
    p_10_temp    = 0;
    p_01_temp    = 0;
    %
    A1_temp  = 0;
    A2_temp  = 0;
    A3_temp  = 0;
    A4_temp  = 0;
    %
    for zz = 1:Sim_times
        % for Markov
        if (snr_SNF_FD(zz) < rho0)
            p_10_temp = p_10_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0)
            p_01_temp = p_01_temp + 1;
        end
        % for F
        % A1
        if (snr_SF_co(zz) >= rho0) && (snr_SF_co(zz) < rho2)
            A1_temp = A1_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) && (snr_SNF_FD(zz) < rho2)
            A2_temp = A2_temp + 1;
        end
        %
        if (snr_SNF_FD(zz) >= rho0) && (snr_SNF_FD(zz) < rho2)
            A3_temp = A3_temp + 1;
        end
        %
        if (snr_SNF_FD(zz) < rho0) && (snr_SF_co(zz) < rho2)
            A4_temp = A4_temp + 1;
        end
    end
    %
    p_10 = p_10_temp/Sim_times;
    p_01 = p_01_temp/Sim_times;
    %
    A1 = A1_temp/Sim_times;
    A2 = A2_temp/Sim_times;
    A3 = A3_temp/Sim_times;
    A4 = A4_temp/Sim_times;
    %
    pi0_sim = p_10/(p_10 + p_01);
    pi1_sim = p_01/(p_10 + p_01);
    %
    %% Analytical results
    a1 = SNR(ss)*pF;
    a2 = SNR(ss)*pN;
    a3 = SNR(ss);
    %
    kappa0 = exp(-rho0/(a1-a2*rho0)/lSN);
    kappa1 = exp(-rho1/a2/lSN);
    kappa2 = exp(-rho2/(a1-a2*rho2)/lSN);
    kappa3 = exp(-max(rho0,rho2)/(a1-a2*max(rho0,rho2))/lSN);
    %
    alpha0 = a3*rho0/(a1-a2*rho0)/lSN;
    alpha1 = a3*rho1/a2/lSN;
    alpha2 = a3*rho2/(a1-a2*rho2)/lSN;
    alpha3 = a3*max(rho0,rho2)/(a1-a2*max(rho0,rho2))/lSN;
    %
    nu1 = (rho1/a2 - rho0/(a1-a2*rho0))*(a1-a2*rho0)/a3/rho0;
    nu2 = (rho2/(a1-a2*rho2) - rho0/(a1-a2*rho0))*(a1-a2*rho0)/a3/rho0;
    %
    mu2 = rho2/(a1-a2*rho2);
    mu1 = rho1/a2;
    %
    % pi_0 and pi_1
    if rho0 < theta
        pi0_ana = Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi)...
            /(Psi(lSF,rho0/(a1-a2*rho0)) + Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi));
        pi1_ana = Psi(lSF,rho0/(a1-a2*rho0))...
            /(Psi(lSF,rho0/(a1-a2*rho0)) + Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi));
    else
        pi0_ana = 1/2;
        pi1_ana = 1/2;
    end
    % F21
    if rho0 >= theta
        F21 = 0;
    else
        if rho2 >= theta
            F21 = Psi(lSF,rho0/(a1-a2*rho0));
        else
            if rho0 < rho2
                F21 = Psi(lSF,rho2/(a1-a2*rho2)) - Psi(lSF,rho0/(a1-a2*rho0));
            else
                F21 = 0;
            end
        end
    end
    % F22
    if rho0 >= theta
        if rho2 < theta
            F22 = Lambda(rho2,kappa2,alpha2,K,a3,lSF,lNF,lrsi);
        else
            F22 = 1;
        end
    else
        if rho2<theta
            F22 = Psi(lSF,rho0/(a1-a2*rho0))...
                - (1-exp(-(rho2/lNF+1/lSF)*rho0/(a1-a2*rho0)))...
                *Omega(rho2,a3,lSF,lNF)*Xi(kappa2,alpha2,K,lrsi);
        else
            F22 = Psi(lSF,rho0/(a1-a2*rho0));
        end
    end
    % F23
    if rho0>theta
        F23 = 0;
    else
        if rho2>=theta
            F23 = 1 - Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi);
        else
            if rho0<rho2
                F23 = Lambda(rho2,kappa2,alpha2,K,a3,lSF,lNF,lrsi)...
                    -Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi);
            else
                F23 = 0;
            end
        end
    end
    % F24
    if rho2 >= theta
        if rho0 < theta
            F24 = Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi);
        else
            F24 = 1;
        end
    else
        if rho0<theta
            F24 = Psi(lSF,rho2/(a1-a2*rho2))...
                - (1-exp(-(rho0/lNF+1/lSF)*mu2))...
                *Omega(rho0,a3,lSF,lNF)*Xi(kappa0,alpha0,K,lrsi);
        else
            F24 = Psi(lSF,rho2/(a1-a2*rho2));
        end
    end
    %
    OP_F_s2_sim(ss) = pi0_sim*(A1+A2) + pi1_sim*(A3+A4);
    OP_F_s2_ana(ss) = pi0_ana*(F21+F22) + pi1_sim*(F23+F24);
end
%% Plot
semilogy(SNR_dB,OP_F_s2_sim,'*')
hold on
semilogy(SNR_dB,OP_F_s2_ana,'-')
%
legend('Simulation','Analysis')
title('OP of User F in Scheme II')
xlabel('SNR (dB)')
ylabel('Outage Probability')
%
save data_OP_F_s2_sim.dat OP_F_s2_sim -ascii
save data_OP_F_s2_ana.dat OP_F_s2_ana -ascii