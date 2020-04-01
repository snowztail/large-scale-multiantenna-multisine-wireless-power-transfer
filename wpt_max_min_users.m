%% * Initialize script for Figure 8
clear; close all; clc; setup; config_max_min_users;

%% * Waveform design by Max-Min-RR and Max-Min-Rand algorithms
minVoltageRr = zeros(length(Variable.nUsers), nRealizations);
minVoltageRand = zeros(length(Variable.nUsers), nRealizations, length(Variable.nCandidates));
for iUser = 1 : length(Variable.nUsers)
    nUsers = Variable.nUsers(iUser);
    weight = ones(1, nUsers);
    [pathloss] = large_scale_fading(distance) * ones(1, nUsers);
    for iRealization = 1 : nRealizations
        channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
        [~, ~, ~, minVoltageRr(iUser, iRealization)] = waveform_max_min_rr(beta2, beta4, txPower, channel, tolerance);
        for iCandidate = 1 : length(Variable.nCandidates)
            nCandidates = Variable.nCandidates(iCandidate);
            [~, ~, ~, minVoltageRand(iUser, iRealization, iCandidate)] = waveform_max_min_rand(beta2, beta4, txPower, channel, tolerance, weight, nCandidates);
        end
    end
end
minVoltageRr = mean(minVoltageRr, 2);
minVoltageRand = squeeze(mean(minVoltageRand, 2));
save('data/wpt_max_min_users.mat');

%% * Result
figure('name', sprintf('Average minimum output voltage of Max-Min-RR and Max-Min-Rand with M = %d and N = %d', nTxs, nSubbands));
for iUser = 1: length(Variable.nUsers)
    subplot(length(Variable.nUsers), 1, iUser);
    stem(1, 1e3 * minVoltageRr(iUser), 'b');
    hold on;
    stem(2 : 1 + length(Variable.nCandidates), 1e3 * minVoltageRand(iUser, :), 'r');
    hold off;
    grid minor;
    legend('Max-Min-RR', 'Max-Min-Rand', 'location', 'nw');
    xticks(1 : 1 + length(Variable.nCandidates));
    xticklabels(string([0 Variable.nCandidates]));
    ylim(1e3 * [min([minVoltageRr(iUser), minVoltageRand(iUser, :)]), max([minVoltageRr(iUser), minVoltageRand(iUser, :)])])
    xlabel('Number of random feasible solutions T');
    ylabel('Average minimum v_{out} [mV]');
    title(sprintf('K = %d', Variable.nUsers(iUser)));
end
savefig('results/wpt_max_min_users.fig');
