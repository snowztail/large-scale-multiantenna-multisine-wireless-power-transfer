clear; close all; clc; initialize; config_wsums;
%% Waveform design by WSum and WSum-S algorithms
voltageWsum = zeros(length(Variable.nUsers), length(Variable.nSubbands), nRealizations);
voltageWSums = zeros(length(Variable.nUsers), length(Variable.nSubbands), nRealizations);
for iUser = 1 : length(Variable.nUsers)
    nUsers = Variable.nUsers(iUser);
    weight = ones(1, nUsers);
    pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10)) * ones(1, nUsers);
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, voltageWsum(iUser, iSubband, iRealization)] = waveform_wsum(beta2, beta4, powerBudget, channel, tolerance, weight);
            [~, voltageWSums(iUser, iSubband, iRealization)] = waveform_wsums(beta2, beta4, powerBudget, channel, tolerance, weight);
        end
    end
end
voltageWsum = mean(voltageWsum, 3);
voltageWSums = mean(voltageWSums, 3);
% save('data/wpt_wsums.mat');
%% Result
legendString = cell(2, length(Variable.nUsers));
legendColor = num2cell(get(gca, 'colororder'), 2);
figure('Name', sprintf('Average output voltage as a function of N with M = 4'));
for iUser = 1 : length(Variable.nUsers)
    plot(Variable.nSubbands, voltageWsum(iUser, :) * 1e3, 'color', legendColor{iUser}, 'Marker', 'o');
    legendString{1, iUser} = sprintf('WSum: K = %d', Variable.nUsers(iUser));
    hold on;
    plot(Variable.nSubbands, voltageWSums(iUser, :) * 1e3, 'color', legendColor{iUser}, 'Marker', 'x');
    legendString{2, iUser} = sprintf('WSum-S: K = %d', Variable.nUsers(iUser));
    hold on;
end
hold off;
xlim([min(Variable.nSubbands), max(Variable.nSubbands)]);
xticks(Variable.nSubbands);
grid on;
legend(legendString(:), 'location', 'nw');
xlabel('Number of tones')
ylabel('Average v_{out} [mV]')
% savefig('results/wpt_wsums.fig');
