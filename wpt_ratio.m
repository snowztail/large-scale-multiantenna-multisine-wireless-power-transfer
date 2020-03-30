clear; close all; clc; initialize; config_ratio;
% Waveform design by CHE Max-Min-Rand and MU UP algorithms
minVoltageCheRand = zeros(length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
minVoltageUp = zeros(length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
for iSubband = 1 : length(Variable.nSubbands)
    nSubbands = Variable.nSubbands(iSubband);
    carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
    for iUser = 1 : length(Variable.nUsers)
        nUsers = Variable.nUsers(iUser);
        weight = ones(1, nUsers);
        pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10)) * ones(1, nUsers);
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, ~, ~, minVoltageCheRand(iSubband, iUser, iRealization)] = waveform_max_min_che_rand(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
            [~, ~, ~, minVoltageUp(iSubband, iUser, iRealization)] = waveform_up(beta2, beta4, powerBudget, channel);
        end
    end
end
minVoltageCheRand = mean(minVoltageCheRand, 3);
minVoltageUp = mean(minVoltageUp, 3);
minVoltageRatio = minVoltageCheRand ./ minVoltageUp;
save('data/wpt_ratio.mat');
%% Result
legendString = cell(length(Variable.nSubbands), 1);
figure('Name', sprintf('Average minimum output voltage ratio: CHE Max-Min-Rand over MU UP with M = %d', nTxs));
hold on;
for iSubband = 1 : length(Variable.nSubbands)
    plot(Variable.nUsers, minVoltageRatio(iSubband, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
end
hold off;
grid minor;
legend(legendString);
xlabel('Number of users');
ylabel('Average minimum output voltage gain ratio');
savefig('results/wpt_ratio.fig');
