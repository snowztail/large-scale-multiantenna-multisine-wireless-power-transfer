%% * Initialize script for Figure 5
clear; close all; clc; setup; config_su_comparison;

%% * Waveform design by SU WPT, CHE WSum, UP and ASS algorithms
voltageSu = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageCheWsum = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageUp = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageAss = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    txPower = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        [carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, voltageSu(iTx, iSubband, iRealization)] = waveform_su(beta2, beta4, txPower, channel, tolerance);
            [~, voltageCheWsum(iTx, iSubband, iRealization)] = waveform_che_wsum(beta2, beta4, txPower, channel, tolerance, weight, pathloss);
            [~, voltageUp(iTx, iSubband, iRealization)] = waveform_up(beta2, beta4, txPower, channel);
            [~, voltageAss(iTx, iSubband, iRealization)] = waveform_ass(beta2, beta4, txPower, channel);
        end
    end
end
voltageSu = mean(voltageSu, 3);
voltageCheWsum = mean(voltageCheWsum, 3);
voltageUp = mean(voltageUp, 3);
voltageAss = mean(voltageAss, 3);
save('data/wpt_su_comparison.mat');

%% * Result
figure('name', sprintf('Average output voltage as a function of (M, N) with K = %d', nUsers));
for iTx = 1: length(Variable.nTxs)
    subplot(length(Variable.nTxs), 1, iTx);
    bar(1e3 * [voltageSu(iTx, :); voltageCheWsum(iTx, :); voltageUp(iTx, :); voltageAss(iTx, :)]');
    grid minor;
    label = [repelem(Variable.nTxs(iTx), length(Variable.nSubbands)); Variable.nSubbands];
    set(gca, 'xticklabel', display_coordinate(label));
    legend('SU WPT', 'CHE WSum', 'UP', 'ASS', 'location', 'nw');
    ylabel('Average v_{out} [mV]');
end
savefig('results/wpt_su_comparison.fig');
