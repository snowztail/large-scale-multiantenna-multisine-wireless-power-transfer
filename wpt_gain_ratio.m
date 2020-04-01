%% * Initialize script for Figure 12
clear; close all; clc; setup; config_gain_ratio;

% * Waveform design by CHE Max-Min-Rand and MU UP algorithms
minVoltageCheRand = zeros(length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
minVoltageUp = zeros(length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
for iSubband = 1 : length(Variable.nSubbands)
    nSubbands = Variable.nSubbands(iSubband);
    [carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
    for iUser = 1 : length(Variable.nUsers)
        nUsers = Variable.nUsers(iUser);
        weight = ones(1, nUsers);
        [pathloss] = large_scale_fading(distance) * ones(1, nUsers);
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, ~, ~, minVoltageCheRand(iSubband, iUser, iRealization)] = waveform_max_min_che_rand(beta2, beta4, txPower, channel, tolerance, pathloss);
            [~, ~, ~, minVoltageUp(iSubband, iUser, iRealization)] = waveform_up(beta2, beta4, txPower, channel);
        end
    end
end
minVoltageCheRand = mean(minVoltageCheRand, 3);
minVoltageUp = mean(minVoltageUp, 3);
minVoltageRatio = minVoltageCheRand ./ minVoltageUp;
save('data/wpt_gain_ratio.mat');

%% * Result
legendString = cell(length(Variable.nSubbands), 1);
figure('name', sprintf('Average minimum output voltage ratio: CHE Max-Min-Rand over MU UP with M = %d', nTxs));
for iSubband = 1 : length(Variable.nSubbands)
    plot(Variable.nUsers, minVoltageRatio(iSubband, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Number of users');
ylabel('Average minimum output voltage gain ratio');
savefig('results/wpt_gain_ratio.fig');
