%% * Initialize script
clear; close all; clc; setup; config_wsum;

%% * Waveform design by SU WPT and WSum algorithms
voltageSu = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageWsum = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    txPower = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        [carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, voltageSu(iTx, iSubband, iRealization)] = waveform_su(beta2, beta4, txPower, channel, tolerance);
            [~, voltageWsum(iTx, iSubband, iRealization)] = waveform_wsum(beta2, beta4, txPower, channel, tolerance, weight);
        end
    end
end
voltageSu = mean(voltageSu, 3);
voltageWsum = mean(voltageWsum, 3);
save('data/wpt_wsum.mat');

%% * Result
figure('name', sprintf('Average output voltage by SU WPT and WSum as a function of sinewaves for single user transmission'));
legendString = cell(2, length(Variable.nTxs));
legendColor = num2cell(get(gca, 'colororder'), 2);
hold on;
for iTx = 1 : length(Variable.nTxs)
    plot(Variable.nSubbands, voltageSu(iTx, :) * 1e3, 'color', legendColor{iTx}, 'marker', 'x');
    legendString{1, iTx} = sprintf('SU WPT: M = %d', Variable.nTxs(iTx));
    plot(Variable.nSubbands, voltageWsum(iTx, :) * 1e3, 'color', legendColor{iTx}, 'marker', 'o');
    legendString{2, iTx} = sprintf('WSum: M = %d', Variable.nTxs(iTx));
end
hold off;
xlim([min(Variable.nSubbands), max(Variable.nSubbands)]);
xticks(Variable.nSubbands);
grid minor;
legend(vec(legendString), 'location', 'nw');
xlabel('Number of tones');
ylabel('Average v_{out} [mV]');
savefig('results/wpt_wsum.fig');
