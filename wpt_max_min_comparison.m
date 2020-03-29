clear; close all; clc; initialize; config_max_min_comparison;
%% Waveform design by Max-Min-Rand, CHE Max-Min-RR and CHE Max-Min-Rand algorithms
minVoltageRand = zeros(length(Variable.nTxs), length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
minVoltageCheRr = zeros(length(Variable.nTxs), length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
minVoltageCheRand = zeros(length(Variable.nTxs), length(Variable.nSubbands), length(Variable.nUsers), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    powerBudget = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
        for iUser = 1 : length(Variable.nUsers)
            nUsers = Variable.nUsers(iUser);
            weight = ones(1, nUsers);
            pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10)) * ones(1, nUsers);
            for iRealization = 1 : nRealizations
                channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
                [~, ~, ~, minVoltageRand(iTx, iSubband, iUser, iRealization)] = waveform_max_min_rand(beta2, beta4, powerBudget, channel, tolerance, weight, nCandidates);
                [~, ~, ~, minVoltageCheRr(iTx, iSubband, iUser, iRealization)] = waveform_max_min_che_rr(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
                [~, ~, ~, minVoltageCheRand(iTx, iSubband, iUser, iRealization)] = waveform_max_min_che_rand(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
            end
        end
    end
end
minVoltageRand = mean(minVoltageRand, 4);
minVoltageCheRr = mean(minVoltageCheRr, 4);
minVoltageCheRand = mean(minVoltageCheRand, 4);
save('data/wpt_max_min_comparison.mat');
%% Result
figure('Name', sprintf('Average minimum output voltage as a function of (M, N, K)'));
bar(1e3 * [vec(minVoltageRand), vec(minVoltageCheRr), vec(minVoltageCheRand)]);
label = [repelem(Variable.nTxs, length(Variable.nSubbands) * length(Variable.nUsers)); repmat(repelem(Variable.nSubbands, length(Variable.nUsers)), [1, length(Variable.nTxs)]); repmat(repmat(Variable.nUsers, [1, length(Variable.nSubbands)]), [1, length(Variable.nTxs)])];
set(gca, 'XTick', 1 : length(Variable.nTxs) * length(Variable.nSubbands) * length(Variable.nUsers), 'xticklabel', display_coordinate(label));
grid on;
legend('Max-Min-Rand', 'CHE Max-Min-RR', 'CHE Max-Min-Rand', 'location', 'nw');
ylabel('Average v_{out} [mV]')
savefig('results/wpt_max_min_comparison.fig');
