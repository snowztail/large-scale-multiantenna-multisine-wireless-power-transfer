clear; close all; clc; initialize; config_su_comparison;
%% Waveform design by SU WPT, CHE WSum, UP and ASS algorithms
voltageSu = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageCheWsum = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageUp = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageAss = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    powerBudget = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nSubbands, nTxs, carrierFrequency, fadingType);
            [~, voltageSu(iTx, iSubband, iRealization)] = waveform_su(beta2, beta4, powerBudget, channel, tolerance);
            [~, voltageCheWsum(iTx, iSubband, iRealization)] = waveform_che_wsum(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
            [~, voltageUp(iTx, iSubband, iRealization)] = waveform_up(beta2, beta4, powerBudget, channel);
            [~, voltageAss(iTx, iSubband, iRealization)] = waveform_ass(beta2, beta4, powerBudget, channel);
        end
    end
end
voltageSu = mean(voltageSu, 3);
voltageCheWsum = mean(voltageCheWsum, 3);
voltageUp = mean(voltageUp, 3);
voltageAss = mean(voltageAss, 3);
save('data/wpt_su_comparison.mat');
%% Result
figure('Name', sprintf('Average output voltage as a function of (M, N) with K = 1'));
for iTx = 1: length(Variable.nTxs)
    subplot(length(Variable.nTxs), 1, iTx);
    bar(1e3 * [voltageSu(iTx, :); voltageCheWsum(iTx, :); voltageUp(iTx, :); voltageAss(iTx, :)]');
    set(gca, 'xticklabel', num2cell(Variable.nSubbands));
    grid on;
    legend('SU WPT', 'CHE WSum', 'UP', 'ASS', 'location', 'nw');
    xlabel('Number of tones');
    ylabel('Average v_{out} [mV]')
    title(sprintf('M = %d', Variable.nTxs(iTx)));
end
savefig('results/wpt_su_comparison.fig');
