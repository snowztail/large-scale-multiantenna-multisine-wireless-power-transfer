clear; close all; clc; setup; config_che_wsum;
%% Waveform design by WSum and CHE WSum algorithms for two-user scenario
userVoltageWsum = zeros(length(Variable.weight), nUsers);
userVoltageCheWsum = zeros(length(Variable.weight), nUsers);
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
for iWeight = 1 : length(Variable.weight)
    weight = Variable.weight(iWeight, :);
    [~, ~, userVoltageWsum(iWeight, :)] = waveform_wsum(beta2, beta4, txPower, channel, tolerance, weight);
    [~, ~, userVoltageCheWsum(iWeight, :)] = waveform_che_wsum(beta2, beta4, txPower, channel, tolerance, weight, pathloss);
end
save('data/wpt_che_wsum.mat');
%% Result
figure('Name', sprintf('Achievable output voltage region with M = 20 and N = 10'));
plot(1e3 * userVoltageWsum(:, 1), 1e3 * userVoltageWsum(:, 2), 'b-');
hold on;
scatter(1e3 * userVoltageCheWsum(:, 1), 1e3 * userVoltageCheWsum(:, 2), 'ro');
hold on;
plot(1e3 * [userVoltageWsum(1, 1), userVoltageWsum(end, 1)], 1e3 * [userVoltageWsum(1, 2), userVoltageWsum(end, 2)], 'b--');
hold on;
plot(1e3 * [userVoltageCheWsum(1, 1), userVoltageCheWsum(end, 1)], 1e3 * [userVoltageCheWsum(1, 2), userVoltageCheWsum(end, 2)], 'r--');
hold off;
grid on;
legend('WSum', 'CHE WSum', 'TDMA: WSum', 'TDMA: CHE WSum', 'location', 'sw');
xlabel('Average v_{out} of user 1 [mV]');
ylabel('Average v_{out} of user 2 [mV]');
savefig('results/wpt_che_wsum.fig');
