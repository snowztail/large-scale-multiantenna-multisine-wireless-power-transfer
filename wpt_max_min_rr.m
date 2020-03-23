clear; close all; clc; initialize; config_max_min_rr;
%% Waveform design by Max-Min-RR and Max-Min-Rand algorithms
minVoltageRr = zeros(length(Variable.nUsers), nRealizations);
minVoltageRand = zeros(length(Variable.nUsers), nRealizations, length(Variable.nCandidates));
for iUser = 1 : length(Variable.nUsers)
    nUsers = Variable.nUsers(iUser);
    weight = ones(1, nUsers);
    pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10)) * ones(1, nUsers);
    for iRealization = 1 : nRealizations
        channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
        [~, ~, ~, minVoltageRr(iUser, iRealization)] = waveform_max_min_rr(beta2, beta4, powerBudget, channel, tolerance, weight);
        for iCandidate = 1 : length(Variable.nCandidates)
            nCandidates = Variable.nCandidates(iCandidate);
            [~, ~, ~, minVoltageRand(iUser, iRealization, iCandidate)] = waveform_max_min_rand(beta2, beta4, powerBudget, channel, tolerance, weight, nCandidates);
        end
    end
end
minVoltageRr = mean(minVoltageRr, 2);
minVoltageRand = squeeze(mean(minVoltageRand, 2));
save('data/wpt_max_min_rr.mat');
%% Result
figure('Name', sprintf('Average minimum output voltage of Max-Min-RR and Max-Min-Rand with M = 4 and R = 4 for K = 2 and 3'));
for iUser = 1: length(Variable.nUsers)
    subplot(length(Variable.nUsers), 1, iUser);
    stem(1e3 * minVoltageRr, 'b');
    hold on;
    stem(1e3 * minVoltageRand, 'r');
    legend('Max-Min-RR', 'Max-Min-Rand');
    grid on;
    ylabel('Average minimum v_{out} [mV]')
    title(sprintf('K = %d', Variable.nUsers(iUser)));
end
savefig('results/wpt_su_comparison.fig');
