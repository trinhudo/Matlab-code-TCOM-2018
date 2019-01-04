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
% Rician parameters
K_dB        = 25; % Rician K-factor in dB
K           = 10.^(K_dB./10);
lrsi        = K/(K+1); % lambda_rsi
% Rician Distribution parameters
chi          = sqrt(K/(K+1)*lrsi); % Noncentrality parameter
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
    B1_temp  = 0;
    B2_temp  = 0;
    B3_temp  = 0;
    B4_temp  = 0;
    B5_temp  = 0;
    B6_temp  = 0;
    B7_temp  = 0;
    B8_temp  = 0;
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
        % for N
        % B1
        if (snr_SNF_FD(zz) >= rho0) && (snr_SN_xF_FD(zz) < rho2)
            B1_temp = B1_temp + 1;
        end
        %
        if (snr_SNF_FD(zz) >= rho0) && (snr_SN_xF_FD(zz) >= rho2) && (snr_SN_xN_FD(zz) < rho1)
            B2_temp = B2_temp + 1;
        end
        %
        if (snr_SNF_FD(zz) < rho0) && (snr_SN_xF_co(zz) < rho2)
            B3_temp = B3_temp + 1;
        end
        %
        if (snr_SNF_FD(zz) < rho0) && (snr_SN_xF_co(zz) >= rho2) && (snr_SN_xN_co(zz) < rho1)
            B4_temp = B4_temp + 1;
        end
        %
        if (snr_SF_co(zz) >= rho0) && (snr_SN_xF_co(zz) < rho2)
            B5_temp = B5_temp + 1;
        end
        %
        if (snr_SF_co(zz) >= rho0) && (snr_SN_xF_co(zz) >= rho2) && (snr_SN_xN_co(zz) < rho1)
            B6_temp = B6_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) && (snr_SN_xF_FD(zz) < rho2)
            B7_temp = B7_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) && (snr_SN_xF_FD(zz) >= rho2) && (snr_SN_xN_FD(zz) <rho1)
            B8_temp = B8_temp + 1;
        end
    end
    %
    p_10 = p_10_temp/Sim_times;
    p_01 = p_01_temp/Sim_times;
    %
    B1 = B1_temp/Sim_times;
    B2 = B2_temp/Sim_times;
    B3 = B3_temp/Sim_times;
    B4 = B4_temp/Sim_times;
    B5 = B5_temp/Sim_times;
    B6 = B6_temp/Sim_times;
    B7 = B7_temp/Sim_times;
    B8 = B8_temp/Sim_times;
    %
    pi0_sim = p_10/(p_10 + p_01);
    pi1_sim = p_01/(p_10 + p_01);
    %
    %     pout_N_s2_sim(ss) = pi_1*(B1 + B2 + B3 + B4) + pi_0*(B5 + B6 + B7 + B8);
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
    mu0 = rho0/(a1-a2*rho0);
    mu2 = rho2/(a1-a2*rho2);
    mu1 = rho1/a2;
    %
    % pi_0 and pi_1
    if rho0 < theta
        pi0_ana = Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi)...
            /(Psi(lSF,mu0) + Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi));
        pi1_ana = Psi(lSF,mu0)...
            /(Psi(lSF,mu0) + Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi));
    else
        pi0_ana = 1/2;
        pi1_ana = 1/2;
    end
    %
    % N21
    if rho0 >= theta
        N21 = 0;
    else
        if rho2<theta
            if (rho0<rho2)
                N21 = (Xi(kappa0,alpha0,K,lrsi) ...
                    - Xi(kappa2,alpha2,K,lrsi))...
                    *Omega(rho0,a3,lSF,lNF);
            else
                N21 = 0;
            end
        else
            N21 = Xi(kappa0,alpha0,K,lrsi)*Omega(rho0,a3,lSF,lNF);
        end
    end
    % N22
    % N22A
    N22A = Omega(rho0,a3,lSF,lNF);
    % N22B
    if (max(rho0,rho2) < theta) && ...
            (max(rho0,rho2)/(a1-a2*max(rho0,rho2)) < rho1/a2)
        N22B = Xi(kappa3,alpha3,K,lrsi) - Xi(kappa1,alpha1,K,lrsi);
    else
        N22B = 0;
    end
    %
    N22 = N22A*N22B;
    % N23
    if rho2 >= theta
        if rho0 < theta
            N23 = Lambda(rho0,kappa0,alpha0,K,a3,lSF,lNF,lrsi);
        else
            N23 = 1;
        end
    else
        N23A = Psi(lSN, rho2/(a1-a2*rho2));
        N23B = Omega(rho0,a3,lSF,lNF);
        % N23C
        if rho0>=theta
            N23C = 0;
        else
            if rho2/(a1-a2*rho2) > rho0/(a1-a2*rho0)
                N23C = Phi(nu2,kappa0,alpha0,K,lrsi) -...
                    Phi(nu2,kappa2,0,K,lrsi);
            else
                N23C = 0;
            end
        end
        %
        N23 = N23A - N23B*N23C;
    end
    % N24
    if (rho2<theta) && (rho2/(a1-a2*rho2) < rho1/a2)
        N24A = exp(-mu2/lSN) - exp(-mu1/lSN);
        N24B = Omega(rho0,a3,lSF,lNF);
        if (rho0>=theta)
            N24C = 0;
        else
            if rho2/(a1-a2*rho2) < rho0/(a1-a2*rho0)
                if rho0/(a1-a2*rho0) < rho1/a2
                    N24C = Phi(nu1,kappa0,alpha0,K,lrsi) ...
                        - Phi(nu1,kappa1,0,K,lrsi);
                else
                    N24C = 0;
                end
            else
                N24C = Phi(nu2,kappa2,0,K,lrsi)...
                    + Phi(nu1,kappa0,alpha0,K,lrsi) ...
                    - Phi(nu1,kappa1,0,K,lrsi)...
                    - Phi(nu2,kappa0,alpha0,K,lrsi);
            end
        end
        %
        N24 = N24A - N24B*N24C;
    else
        N24 = 0;
    end
    % N25
    % N25A
    if rho0 < theta
        N25A = 1 - Psi(lSF,rho0/(a1-a2*rho0));
    else
        N25A = 0;
    end
    % N25B
    if rho2 < theta
        N25B = Psi(lSN,rho2/(a1-a2*rho2));
    else
        N25B = 1;
    end
    %
    N25 = N25A*N25B;
    % N26
    % N26A = N25A
    if rho0 < theta
        N26A = 1 - Psi(lSF,rho0/(a1-a2*rho0));
    else
        N26A = 0;
    end
    % N26B
    if (rho2 < theta) && (rho2/(a1-a2*rho2) < rho1/a2)
        N26B = Psi(lSN,rho1/a2) - Psi(lSN,rho2/(a1-a2*rho2));
    else
        N26B = 0;
    end
    N26 = N26A*N26B;
    % N27
    % N27A
    if rho0 < theta
        N27A = Psi(lSF,rho0/(a1-a2*rho0));
    else
        N27A = 1;
    end
    % N27B
    if rho2 < theta
        N27B = 1 - Xi(kappa2,alpha2,K,lrsi);
    else
        N27B = 1;
    end
    N27 = N27A*N27B;
    % N28
    % N28A
    if rho0 < theta
        N28A = Psi(lSF,rho0/(a1-a2*rho0));
    else
        N28A = 1;
    end
    % N28B
    if (rho2<theta) && (rho1/a2 > rho2/(a1-a2*rho2))
        N28B = Xi(kappa2,alpha2,K,lrsi) - Xi(kappa1,alpha1,K,lrsi);
    else
        N28B = 0;
    end
    N28 = N28A*N28B;
    %
    OP_N_s2_sim(ss) = pi1_sim*(B1 + B2 + B3 + B4) ...
        + pi0_sim*(B5 + B6 + B7 + B8);
    OP_N_s2_ana(ss) = pi1_ana*(N21+ N22+ N23 + N24) ...
        + pi0_ana*(N25+N26+N27+N28);
    %
    %     test_sim(ss) = B4;
    %     test_ana(ss) = N24;
end
%% Plot
% semilogy(SNR_dB,test_sim,'o',SNR_dB,test_ana,'-')
semilogy(SNR_dB,OP_N_s2_sim,'o',SNR_dB,OP_N_s2_ana,'-')
%
legend('Simulation','Analysis')
title('OP of User N in Scheme II')
xlabel('SNR (dB)')
ylabel('Outage Probability')
%
save data_OP_N_s2_sim.dat OP_N_s2_sim -ascii
save data_OP_N_s2_ana.dat OP_N_s2_ana -ascii