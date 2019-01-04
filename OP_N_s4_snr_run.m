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
R_on_off    = 1; % for optimal ON/OFF
rho2        = 2^RtF - 1; % full time slot
rho1        = 2^RtN - 1;
rho2h       = 2^(2*RtF) - 1; % half time slot
rho1h       = 2^(2*RtN) - 1;
rho0        = 2.^R_on_off - 1; % SNR threshold for ON/OFF
rho0h        = 2.^(2*R_on_off) - 1; % SNR threshold for ON/OFF

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
        if (snr_SNF_HD(zz) < rho0h)
            p_10_temp = p_10_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0)
            p_01_temp = p_01_temp + 1;
        end
        % for N
        % B1
        if (snr_SNF_HD(zz) >= rho0h) && (snr_SN_xF_co(zz) < rho2h)
            B1_temp = B1_temp + 1;
        end
        % B2
        if (snr_SNF_HD(zz) >= rho0h) && (snr_SN_xF_co(zz) >= rho2h) ...
                && (snr_SN_xN_co(zz) < rho1h)
            B2_temp = B2_temp + 1;
        end
        %
        if (snr_SNF_HD(zz) < rho0h) && (snr_SN_xF_co(zz) < rho2)
            B3_temp = B3_temp + 1;
        end
        %
        if (snr_SNF_HD(zz) < rho0h) && (snr_SN_xF_co(zz) >= rho2) ...
                && (snr_SN_xN_co(zz) < rho1)
            B4_temp = B4_temp + 1;
        end
        %
        if (snr_SF_co(zz) >= rho0) && (snr_SN_xF_co(zz) < rho2)
            B5_temp = B5_temp + 1;
        end
        %
        if (snr_SF_co(zz) >= rho0) && (snr_SN_xF_co(zz) >= rho2) ...
                && (snr_SN_xN_co(zz) < rho1)
            B6_temp = B6_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) && (snr_SN_xF_co(zz) < rho2h)
            B7_temp = B7_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) && (snr_SN_xF_co(zz) >= rho2h) ...
                && (snr_SN_xN_co(zz) <rho1h)
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
    pi_0_sim = p_10/(p_10 + p_01);
    pi_1_sim = p_01/(p_10 + p_01);
    %
    %% Analysis
    a1 = SNR(ss)*pF;
    a2 = SNR(ss)*pN;
    a3 = SNR(ss);
    %
    chi = 1/(a2*a3*lNF*lSF);
    xi = a1/a2;
    mu0     = rho0/(a1-a2*rho0);
    mu0h    = rho0h/(a1-a2*rho0h);
    mu1     = rho1/a2;
    mu1h    = rho1h/a2;
    mu2     = rho2/(a1-a2*rho2);
    mu2h    = rho2h/(a1-a2*rho2h);
    mu3h    = max(rho0h,rho2h)/(a1-a2*max(rho0h,rho2h));
    %
    Upsilon_mu0h_rho0h = 1 - exp(-rho0h/lSF/(a1-a2*rho0h))  ...
        - 1/a2/a3/lNF/lSF*exp(-rho0h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*rho0h/(a1-a2*rho0h)+a3*lNF,chi,xi));
    % pi_0 and pi_1
    if (rho0 < theta) && (rho0h < theta)
        pi_0_ana = (1-(1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h))/ ...
            (Psi(lSF,mu0) + (1-(1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h)));
        pi_1_ana = Psi(lSF,mu0)/ ...
            (Psi(lSF,mu0) + (1-(1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h)));
    elseif (rho0 < theta) && (rho0h >= theta)
        pi_0_ana = 1/(1+Psi(lSF,mu0));
        pi_1_ana = Psi(lSF,mu0)/(1+Psi(lSF,mu0));
    elseif (rho0 >= theta) && (rho0h >= theta)
        pi_0_ana = 1/2;
        pi_1_ana = 1/2;
    end
    %
    % N41
    if (rho0h >= theta)
        N41 = 0;
    else
        if (rho2h >= theta)
            N41 = (1-Upsilon_mu0h_rho0h)*(1-Psi(lSN,mu0h));
        else
            if (rho0h < rho2h)
                N41 = (1-Upsilon_mu0h_rho0h)*(Psi(lSN,mu2h)-Psi(lSN,mu0h));
            else
                N41 = 0;
            end
        end
    end
    % N42
    if (rho0h < theta) && (rho2h < theta)
        if mu3h < mu1h
            N42 = (1-Upsilon_mu0h_rho0h)*(Psi(lSN,mu1h) - Psi(lSN,mu3h));
        else
            N42 = 0;
        end
    else
        N42 = 0;
    end
    % N43
    if rho2 >= theta
        if rho0h < theta
            N43 = 1 - (1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h);
        else
            N43 = 1;
        end
    else
        if rho0h >= theta
            N43 = Psi(lSN,mu2);
        else
            if mu0h < mu2
                N43 = Psi(lSN,mu2) ...
                    - (1-Upsilon_mu0h_rho0h)*(exp(-mu0h/lSN)-exp(-mu2/lSN));
            else
                N43 = Psi(lSN,mu2);
            end
        end
    end
    % N44
    if rho2 >= theta
        N44 = 0;
    else
        if mu1 <= mu2
            N44= 0;
        else
            if rho0h >= theta
                N44 = Psi(lSN,mu1)-Psi(lSN,mu2);
            else
                if mu0h < mu2
                    N44 = Psi(lSN,mu1)-Psi(lSN,mu2) ...
                        - (1-Upsilon_mu0h_rho0h)*(Psi(lSN,mu1)-Psi(lSN,mu2));
                elseif mu2 <= mu0h && mu0h<=mu1
                    N44 = Psi(lSN,mu1)-Psi(lSN,mu2) ...
                        - (1-Upsilon_mu0h_rho0h)*(Psi(lSN,mu1)-Psi(lSN,mu0h));
                elseif mu0h>=mu1
                    N44 = Psi(lSN,mu1)-Psi(lSN,mu2);
                end
            end
        end
    end
    % N45
    if rho0>=theta
        N45 = 0;
    else
        if rho2>=theta
            N45 = 1 - Psi(lSF,mu0);
        else
            N45 = (1-Psi(lSF,mu0))*Psi(lSN,mu2);
        end
    end
    % N46
    if (rho0>=theta) || (rho2>=theta)
        N46 = 0;
    elseif (rho0<theta) && (rho2<theta)
        if mu1<=mu2
            N46 = 0;
        else
            N46 = (1-Psi(lSF,mu0))*(Psi(lSN,mu1)-Psi(lSN,mu2));
        end
    end
    % N47
    if rho0<theta && rho2h<theta
        N47 = Psi(lSF,mu0)*Psi(lSN,mu2h);
    elseif rho0<theta && rho2h>=theta
        N47 = Psi(lSF,mu0);
    elseif rho0>=theta && rho2h<theta
        N47 = Psi(lSN,mu2h);
    elseif rho0>theta && rho2h>=theta
        N47 = 1;
    end
    % N48
    if rho2h>=theta
        N48=0;
    else
        if mu1h<=mu2h
            N48 = 0;
        else
            if rho0<theta
                N48 = Psi(lSF,mu0)*(Psi(lSN,mu1h)-Psi(lSN,mu2h));
            else
                N48 = Psi(lSN,mu1h)-Psi(lSN,mu2h);
            end
        end
    end
    OP_N_s4_sim(ss) = pi_1_sim*(B1 + B2 + B3 + B4) ...
        + pi_0_sim*(B5 + B6 + B7 + B8);
    OP_N_s4_ana(ss) = pi_1_ana*(N41+N42+N43+N44) ...
        + pi_0_ana*(N45+N46+N47+N48);
    %
    %     test_sim = B4;
    %     test_ana = N42;
end
%% Plot
% semilogy(SNR_dB,test_sim,'bo:',SNR_dB,test_ana,'b-')
semilogy(SNR_dB,OP_N_s4_sim,'bo:',SNR_dB,OP_N_s4_ana,'b-')
%
legend('Simulation','Analysis')
%
title('SCHEME IV, USER N')
xlabel('SNR (dB)')
ylabel('Outage Probability')
%
save data_OP_N_s4_sim.dat OP_N_s4_sim -ascii
save data_OP_N_s4_ana.dat OP_N_s4_ana -ascii