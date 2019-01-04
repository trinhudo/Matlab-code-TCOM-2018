close all
clear all
%
SNR_dB = -10:5:20; % in dB
%
load data_OP_F_co_sim.dat
load data_OP_F_co_ana.dat
load data_OP_F_s1_sim.dat
load data_OP_F_s1_ana.dat
load data_OP_F_s2_sim.dat
load data_OP_F_s2_ana.dat
load data_OP_F_s3_sim.dat
load data_OP_F_s3_ana.dat
load data_OP_F_s4_sim.dat
load data_OP_F_s4_ana.dat
%
semilogy(SNR_dB,data_OP_F_co_sim,'ko');
hold on
semilogy(SNR_dB,data_OP_F_s1_sim,'b*');
hold on
semilogy(SNR_dB,data_OP_F_s2_sim,'rs');
hold on
semilogy(SNR_dB,data_OP_F_s3_sim,'gv');
hold on
semilogy(SNR_dB,data_OP_F_s4_sim,'cp');
hold on
semilogy(SNR_dB,data_OP_F_co_ana,'k-');
hold on
semilogy(SNR_dB,data_OP_F_s1_ana,'b-');
hold on
semilogy(SNR_dB,data_OP_F_s2_ana,'r-');
hold on
semilogy(SNR_dB,data_OP_F_s3_ana,'g-');
hold on
semilogy(SNR_dB,data_OP_F_s4_ana,'c-');
hold off
axis([-10 20 1e-4 1e0])
%
legend('Non-coop. NOMA (sim.)', ...
    'only-FDR (sim.)', ...
    'on/off-FDR (sim.)', ...
    'only-HDR (sim.)', ...
    'on/off-HDR (sim.)', ...
    'Location','SouthWest')
%
xlabel('Transmit SNR (dBm)')
ylabel('Outage Probability (OP)')
%
xticks([-10 0 10 20])
xticklabels([20 30 40 50])
% set(gca, 'LooseInset', [0,0,0,0]);