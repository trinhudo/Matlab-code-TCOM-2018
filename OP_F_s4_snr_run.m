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
    A1_temp  = 0;
    A2_temp  = 0;
    A3_temp  = 0;
    A4_temp  = 0;
    A5_temp  = 0;
    A6_temp  = 0;
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
        % for User F
        % A1
        if (snr_SF_co(zz) >= rho0) && (snr_SF_co(zz) < rho2)
            A1_temp = A1_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) &&...
                (snr_SN_xF_co(zz)<rho2h) && (snr_SF_co(zz)<rho2h)
            A2_temp = A2_temp + 1;
        end
        %
        if (snr_SF_co(zz) < rho0) &&...
                (snr_SN_xF_co(zz)>=rho2h) && ...
                (snr_SF_co(zz)+snr_NF_HD(zz)<rho2h)
            A3_temp = A3_temp + 1;
        end
        %
        if (snr_SNF_HD(zz) >= rho0h) &&...
                (snr_SN_xF_co(zz)<rho2h) && (snr_SF_co(zz)<rho2h)
            A4_temp = A4_temp + 1;
        end
        %
        if (snr_SNF_HD(zz) >= rho0h) &&...
                (snr_SN_xF_co(zz)>=rho2h) && (snr_SF_co(zz)+snr_NF_HD(zz)<rho2h)
            A5_temp = A5_temp + 1;
        end
        %
        if (snr_SNF_HD(zz) < rho0h) && (snr_SF_co(zz) < rho2)
            A6_temp = A6_temp + 1;
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
    A5 = A5_temp/Sim_times;
    A6 = A6_temp/Sim_times;
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
    mu0h     = rho0h/(a1-a2*rho0h);
    mu1     = rho1/a2;
    mu1h    = rho1h/a2;
    mu2     = rho2/(a1-a2*rho2);
    mu2h    = rho2h/(a1-a2*rho2h);
    mu3h    = max(rho0,rho2h)/(a1-a2*max(rho0,rho2h));
    rho3h   = min(rho0,rho2h);
    mu02    = min(mu0,mu2);
    mu02hat = min(mu0,mu2h);
    mu02tilde = min(mu0h,mu2h);
    mu02check = min(mu0h,mu2);
    %
    Upsilon_mu0h_rho0h = 1 - exp(-rho0h/lSF/(a1-a2*rho0h))  ...
        - 1/a2/a3/lNF/lSF*exp(-rho0h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*rho0h/(a1-a2*rho0h)+a3*lNF,chi,xi));
    %
    Upsilon_mu2h_rho2h = 1 - exp(-mu2h/lSF)  ...
        - 1/a2/a3/lNF/lSF*exp(-rho2h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*mu2h+a3*lNF,chi,xi));
    %
    Upsilon_mu02hat_rho2h = 1 - exp(-mu02hat/lSF)  ...
        - 1/a2/a3/lNF/lSF*exp(-rho2h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*mu02hat+a3*lNF,chi,xi));
    %
    Upsilon_mu02tilde_rho0h = 1 - exp(-mu02tilde/lSF)  ...
        - 1/a2/a3/lNF/lSF*exp(-rho0h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*mu02tilde+a3*lNF,chi,xi));
    %
    Upsilon_mu02check_rho0h = 1 - exp(-mu02check/lSF)  ...
        - 1/a2/a3/lNF/lSF*exp(-rho0h/a3/lNF + a1/a2/a3/lNF + 1/a2/lSF) ...
        *(ApproxIntegral(a3*lNF,chi,xi) ...
        - ApproxIntegral(a2*a3*lNF*mu02check+a3*lNF,chi,xi));
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
    % F41
    if rho0 >= theta
        F41 = 0;
    else
        if rho2 >= theta
            F41 = 1 - Psi(lSF,mu0);
        else
            if rho0 < rho2
                F41 = Psi(lSF,mu2) - Psi(lSF,mu0);
            else
                F41 = 0;
            end
        end
    end
    % F42
    if rho0 >= theta && rho2h < theta
        F42 = Psi(lSF,mu2h)*Psi(lSN,mu2h);
    elseif rho0 >= theta && rho2h >= theta
        F42 = 1;
    elseif rho0 < theta && rho2h >= theta
        F42 = Psi(lSF,mu0);
    elseif rho0 < theta && rho2h < theta
        F42 = Psi(lSF,mu02hat)*Psi(lSN,mu2h);
    end
    % F43
    if rho2h >= theta
        F43 = 0;
    else
        if rho0 >= theta
            F43 = (1-Psi(lSN,mu2h))*Upsilon_mu2h_rho2h;
        else
            F43 = (1-Psi(lSN,mu2h))*Upsilon_mu02hat_rho2h;
        end
    end
    % F44
    if rho2h >= theta
        F44 = 1 - (1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h);
    else
        if rho0h >= theta
            F44 = 0;
        else
            if rho0h >= rho2h
                F44 = 0;
            else
                F44 = (Psi(lSN,mu2h)-Psi(lSN,mu0h)) ...
                    *(Psi(lSF,mu2h) - Upsilon_mu02tilde_rho0h);
            end
        end
    end
    % F45
    if (rho0h >= theta) || (rho2h >= theta)
        F45 = 0;
    else
        if rho0h >= rho2h
            F45 = 0;
        else
            F45 = (1-Psi(lSN,mu2h)) ...
                * (Upsilon_mu2h_rho2h - Upsilon_mu0h_rho0h);
        end
    end
    % F46
    if rho2 >= theta
        if rho0h >= theta
            F46 = 1;
        else
            F46 = 1 - (1-Psi(lSN,mu0h))*(1-Upsilon_mu0h_rho0h);
        end
    else
        if rho0h >= theta
            F46 = Psi(lSF,mu2);
        else
            F46 = Psi(lSF,mu2) ...
                - (1-Psi(lSN,mu0h))*(Psi(lSF,mu2)-Upsilon_mu02check_rho0h);
        end
    end
    %
    OP_F_s4_sim(ss) = pi_0_sim*(A1 + A2 + A3) + pi_1_sim*(A4 + A5 + A6);
    OP_F_s4_ana(ss) = pi_0_ana*(F41+F42+F43) + pi_1_ana*(F44+F45+F46);
    %%
end
%% Plot
semilogy(SNR_dB,OP_F_s4_sim,'rs:')
hold on
semilogy(SNR_dB,OP_F_s4_ana,'r-');
%
legend('Simulation','Analysis')
title('SCHEME IV, USER F')
%
xlabel('SNR (dB)')
ylabel('Outage Probability')
%
save data_OP_F_s4_sim.dat OP_F_s4_sim -ascii
save data_OP_F_s4_ana.dat OP_F_s4_ana -ascii