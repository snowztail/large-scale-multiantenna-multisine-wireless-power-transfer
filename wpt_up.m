clear; close all; clc; setup; config_up;
%% Waveform design by CHE Max-Min-Rand and MU UP algorithms
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
%             [~, ~, ~, minVoltageCheRand(iSubband, iUser, iRealization)] = waveform_max_min_che_rand(beta2, beta4, txPower, channel, tolerance, weight, pathloss);
            [~, ~, ~, minVoltageUp(iSubband, iUser, iRealization)] = waveform_up(beta2, beta4, txPower, channel);
        end
    end
end
minVoltageCheRand = mean(minVoltageCheRand, 3);
minVoltageUp = mean(minVoltageUp, 3);
save('data/wpt_up.mat');
%% Result
figure('Name', sprintf('Average minimum output voltage as a function of (N, K) with M = %d', nTxs));
bar(1e3 * [vec(minVoltageCheRand), vec(minVoltageUp)]);
label = [repelem(Variable.nSubbands, length(Variable.nUsers)); repmat(Variable.nUsers, [1, length(Variable.nSubbands)])];
set(gca, 'xticklabel', display_coordinate(label));
grid on;
legend('CHE Max-Min-Rand', 'MU UP', 'location', 'ne');
ylabel('Average minimum v_{out} [mV]')
savefig('results/wpt_up.fig');
