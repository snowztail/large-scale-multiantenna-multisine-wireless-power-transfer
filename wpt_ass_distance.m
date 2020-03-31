%% * Initialize script for Figure 3b
clear; close all; clc; setup; config_ass_distance;

%% * Waveform design by ASS algorithm
voltage = zeros(length(Variable.nTxs), length(Variable.distance), nRealizations);
for iCase = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iCase);
    nSubbands = Variable.nSubbands(iCase);
    [carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
    for iDistance = 1 : length(Variable.distance)
        distance = Variable.distance(iDistance);
        [pathloss] = large_scale_fading(distance) * ones(nUsers, 1);
        for iRealization = 1 : nRealizations
            channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
            [~, voltage(iCase, iDistance, iRealization)] = waveform_ass(beta2, beta4, txPower, channel);
        end
    end
end
voltage = mean(voltage, 3);
save('data/wpt_ass_distance.mat');

%% * Result
figure('name', sprintf('Average output voltage by ASS as a function of distance'));
legendString = cell(length(Variable.nTxs), 1);
for iCase = 1 : length(Variable.nTxs)
    semilogy(Variable.distance, voltage(iCase, :));
    legendString{iCase} = sprintf('M = %d, N = %d', Variable.nTxs(iCase), Variable.nSubbands(iCase));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Distance [m]');
ylabel('Average v_{out} [V]');
savefig('results/wpt_ass_distance.fig');
