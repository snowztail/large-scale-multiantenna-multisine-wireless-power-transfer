clear; close all; clc; setup; config_cdf;
%% Waveform design by EQ WSum, FA WSum, Max-Min-Rand, and CHE Max-Min-Rand algorithms
voltageWsumEq = cell(length(Variable.nUsers), nRealizations);
voltageWsumFa = cell(length(Variable.nUsers), nRealizations);
voltageRand = cell(length(Variable.nUsers), nRealizations);
voltageCheRand = cell(length(Variable.nUsers), nRealizations);
for iUser = 1 : length(Variable.nUsers)
    nUsers = Variable.nUsers(iUser);
    weightEq = ones(1, nUsers);
    [pathloss] = large_scale_fading(distance) * ones(1, nUsers);
    for iRealization = 1 : nRealizations
        channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
        % assign FA weight based on voltage achieved by UP
        [~, ~, voltageUp] = waveform_up(beta2, beta4, txPower, channel);
        weightFa = voltageUp .^ (-1) / sum(voltageUp .^ (-1));
        [~, ~, voltageWsumEq{iUser, iRealization}] = waveform_wsum(beta2, beta4, txPower, channel, tolerance, weightEq);
        [~, ~, voltageWsumFa{iUser, iRealization}] = waveform_wsum(beta2, beta4, txPower, channel, tolerance, weightFa);
        [~, ~, voltageRand{iUser, iRealization}] = waveform_max_min_rand(beta2, beta4, txPower, channel, tolerance, weightEq, nCandidates);
        [~, ~, voltageCheRand{iUser, iRealization}] = waveform_max_min_che_rand(beta2, beta4, txPower, channel, tolerance, weightEq, pathloss);
    end
end
save('data/wpt_cdf.mat');
%% Result
figure('name', sprintf('CDF of output voltage with M = %d and N = %d', nTxs, nSubbands));
legendString = cell(4, length(Variable.nUsers));
legendColor = num2cell(get(gca, 'colororder'), 2);
for iUser = 1 : length(Variable.nUsers)
    plotWsumEq = cdfplot(1e3 * cell2mat(voltageWsumEq(iUser, :)));
    set(plotWsumEq, 'color', legendColor{iUser}, 'LineStyle', ':');
    legendString{1, iUser} = sprintf('WSum-EQ: K = %d', Variable.nUsers(iUser));
    hold on;
    plotWsumFa = cdfplot(1e3 * cell2mat(voltageWsumFa(iUser, :)));
    set(plotWsumFa, 'color', legendColor{iUser}, 'LineStyle', '-.');
    legendString{2, iUser} = sprintf('WSum-FA: K = %d', Variable.nUsers(iUser));
    hold on;
    plotRand = cdfplot(1e3 * cell2mat(voltageRand(iUser, :)));
    set(plotRand, 'color', legendColor{iUser}, 'LineStyle', '-');
    legendString{3, iUser} = sprintf('Max-Min-Rand: K = %d', Variable.nUsers(iUser));
    hold on;
    plotCheRand = cdfplot(1e3 * cell2mat(voltageCheRand(iUser, :)));
    set(plotCheRand, 'color', legendColor{iUser}, 'LineStyle', '--');
    legendString{4, iUser} = sprintf('CHE Max-Min-Rand: K = %d', Variable.nUsers(iUser));
    hold on;
end
hold off;
grid minor;
set(gca, 'XScale', 'log')
legend(legendString(:), 'location', 'nw');
xlabel('Average v_{out} [mV]');
ylabel('CDF');
savefig('results/wpt_cdf.fig');
