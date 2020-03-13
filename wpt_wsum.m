clear; close all; clc; initialize; config_wsum;
%% Waveform design by SU WPT and WSum algorithms
voltageSu = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageWsum = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    powerBudget = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
        for iRealization = 1 : nRealizations
            % \boldsymbol{h}_{n}
            channel = channel_tgn_e(pathloss, nSubbands, nTxs, carrierFrequency, fadingType);
            % \boldsymbol{s_n}, v_{\text{out},q}
            [~, voltageSu(iTx, iSubband, iRealization)] = waveform_su(beta2, beta4, powerBudget, channel, tolerance);
            [~, voltageWsum(iTx, iSubband, iRealization)] = waveform_wsum(beta2, beta4, powerBudget, channel, tolerance, weight);
        end
    end
end
voltageSu = mean(voltageSu, 3);
voltageWsum = mean(voltageWsum, 3);
save('data/wpt_wsum.mat');
%% Result
legendString = cell(2, length(Variable.nTxs));
legendColor = num2cell(get(gca, 'colororder'), 2);
figure('Name', sprintf('Average single user output voltage by SU WPT and WSum as a function of sinewaves'));
for iTx = 1 : length(Variable.nTxs)
    plot(Variable.nSubbands, voltageSu(iTx, :) * 1e3, 'color', legendColor{iTx}, 'Marker', 'x');
    legendString{1, iTx} = sprintf('SU WPT: M = %d', Variable.nTxs(iTx));
    hold on;
    plot(Variable.nSubbands, voltageWsum(iTx, :) * 1e3, 'color', legendColor{iTx}, 'Marker', 'o');
    legendString{2, iTx} = sprintf('WSum: M = %d', Variable.nTxs(iTx));
    hold on;
end
hold off;
xlim([min(Variable.nSubbands), max(Variable.nSubbands)]);
xticks(Variable.nSubbands);
grid on;
legend(legendString(:), 'location', 'nw');
xlabel('Number of tones')
ylabel('Average v_{out} [mV]')
savefig('results/wpt_wsum.fig');
